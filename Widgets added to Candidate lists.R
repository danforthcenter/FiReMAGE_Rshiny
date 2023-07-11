library(shiny)
library(data.table)
library(doParallel)
library(reader)
library(ComplexUpset)
library(ggplot2)
library(ggrepel)
library(zip)
library(dplyr)


# Define UI for app
ui <- fluidPage(
  titlePanel("Candidate Lists"),
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload CSV Files or ZIP file containing CSVs",
                multiple = TRUE,
                accept = c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv',
                           '.zip')),
      uiOutput("species_select"), 
      uiOutput("trait_select")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Violin Plot", plotOutput("candidatePlot")),
        tabPanel("Upset Plot", plotOutput("upsetPlot")),
        tabPanel("Summary Table", tableOutput("summaryTable"))
      )
    )
  )
)

# Define server logic 
server <- function(input, output) {
  processed_files <- reactive({
    file_paths <- NULL
    for (i in seq_along(input$files$name)) {
      ext <- tools::file_ext(input$files$name[[i]])
      if (ext == "zip") {
        temp_dir <- tempdir()
        unzip(input$files$datapath[[i]], exdir = temp_dir)
        csv_files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE, pattern = "\\.csv$")
        file_paths <- c(file_paths, csv_files)
      } else {
        file_paths <- c(file_paths, input$files$datapath[[i]])
      }
    }
    file_paths
  })
  
  data_list <- reactive({
    lapply(processed_files(), function(fp) {
      if (!dir.exists(fp)) {
        data.table::fread(fp, sep = ",", header = T, select = c(
          "org",
          "loci",
          "Orthogroup",
          "trait",
          "Gene Name",
          "Description",
          "present",
          "genecount",
          "pFDR"
        ))
      }
    })
  })
  
  candidateList <- reactive({
    if (is.null(data_list())) return(NULL)
    rbindlist(data_list(), fill = TRUE)
  })
  
  # Creating the dropdown for species selection
  output$species_select <- renderUI({
    if (is.null(candidateList())) return(NULL)
    candidateList <- candidateList()
    selectInput("selected_species", "Select Species", choices = unique(candidateList$org), selected = unique(candidateList$org)[1], multiple = TRUE)
  })
  
  # Creating the dropdown for trait selection
  output$trait_select <- renderUI({
    if (is.null(candidateList())) return(NULL)
    candidateList <- candidateList()
    selectInput("selected_trait", "Select Trait", choices = unique(candidateList$trait), selected = unique(candidateList$trait)[1], multiple = TRUE)
  })
  
  output$candidatePlot <- renderPlot({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    candidateList <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
    if (nrow(candidateList) == 0) return(NULL)
    candidateList$present <- gsub("^([0-9])$", "\\1/5", candidateList$present)
    candidateList$present <- factor(candidateList$present, levels = c("5/5", "4/5", "3/5"))
    trait_labels <- aggregate(candidateList$pFDR, list(candidateList$trait, candidateList$present), FUN = mean)
    colnames(trait_labels) <- c("trait", "present", "pFDR")
    
    ggplot(candidateList, aes(x = present, y = pFDR, fill = present)) +
      geom_violin(color = NA) +
      geom_text_repel(
        data = trait_labels,
        aes(label = trait),
        size = 4.5,
        nudge_y = 0.08,
        direction = "both",
        min.segment.length = 2,
        force = 3,
        max.overlaps = 40
      ) +
      scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb"), guide = "none") +
      theme_bw(base_size = 18) +
      theme(axis.text.x = element_text(angle = 40, hjust = 1)) +
      ylab("FDR: mean(Permutations)/Actual") + xlab("species in orthogroup")
  })
  
  
  output$upsetPlot <- renderPlot({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    candidateList <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
    casted_groups <- dcast(candidateList, Orthogroup + trait ~ org, fill = 0, value.var = "loci", fun.aggregate = function(x) {
      length(x) > 0
    })
    
    upset(
      casted_groups[, -c(1, 2)],
      as.character(colnames(casted_groups)[-c(1, 2)]),
      name = "Orthogroup species composition",
      sort_intersections_by = c("degree", "cardinality"),
      sort_sets = FALSE,
      stripes = c("white", "white"),
      matrix = (
        intersection_matrix(geom = geom_point(size = 4),
                            segment = geom_segment(linetype = NA))
      ),
      base_annotations = list(
        'Intersection size' = intersection_size(
          colour = "black",
          text_colors = c("red", "blue"),  # Changing color of numbers
          text = list(size = 6)
        )
      ),
      themes = upset_default_themes(text = element_text(size = 25))
    )
  })
  
  output$summaryTable <- renderTable({
    if (is.null(candidateList())) return()
    candidateList <- candidateList()
    summary_stats <- candidateList[, .(Mean = mean(pFDR), SD = sd(pFDR)), by = trait]
    summary_stats
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
