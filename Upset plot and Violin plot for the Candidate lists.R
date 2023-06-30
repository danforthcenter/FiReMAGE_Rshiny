# Load necessary libraries
library(shiny)
library(data.table)
library(doParallel)
library(reader)
library(ComplexUpset)
library(ggplot2)
library(ggrepel)
library(zip)

# UI for the app
ui <- fluidPage(
  titlePanel("Candidate Lists"),
  sidebarLayout(
    sidebarPanel(
      # Allow the user to upload either a single CSV, multiple CSVs or a zip file containing CSVs
      fileInput("files", "Upload CSV Files or ZIP file containing CSVs",
                multiple = TRUE,
                accept = c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv',
                           '.zip'))
    ),
    mainPanel(
      # Main panel consists of two tabs - one for Violin Plot and one for Upset Plot
      tabsetPanel(
        tabPanel("Violin Plot", plotOutput("candidatePlot")),
        tabPanel("Upset Plot", plotOutput("upsetPlot"))
      )
    )
  )
)

# Server function
server <- function(input, output) {
  # Process the uploaded files - extract CSVs from ZIP files if necessary
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
  
  # Read all uploaded data into a list of data.tables
  data_list <- reactive({
    lapply(processed_files(), function(fp) {
      if (!dir.exists(fp)) {
        data.table::fread(fp, sep = ",", header = T)
      }
    })
  })
  
  # Combine all uploaded data into one table
  candidateList <- reactive({
    if (is.null(data_list())) return(NULL)
    rbindlist(data_list(), fill = TRUE)
  })
  
  # Generate the Violin plot
  output$candidatePlot <- renderPlot({
    if (is.null(candidateList())) return()
    candidateList <- candidateList()
    if (nrow(candidateList) == 0) return(NULL)
    
    # Preprocess the data for plotting
    candidateList$present <- gsub("^([0-9])$", "\\1/5", candidateList$present)
    candidateList$present <- factor(candidateList$present, levels = c("5/5", "4/5", "3/5"))
    
    # Compute the average pFDR per trait and presence
    trait_labels <- aggregate(candidateList$pFDR, list(candidateList$trait, candidateList$present), FUN = mean)
    colnames(trait_labels) <- c("trait", "present", "pFDR")
    
    # Generate the Violin plot with text repelling for better visibility of labels
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
  
  # Generate the Upset plot
  output$upsetPlot <- renderPlot({
    if (is.null(candidateList())) return()
    candidateList <- candidateList()
    
    # Cast the data to the required format for Upset plotting
    casted_groups <- dcast(candidateList, Orthogroup + trait ~ org, fill = 0, value.var = "loci", fun.aggregate = function(x) {
      length(x) > 0
    })
    
    # Generate the Upset plot
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
          text_colors = c("black", "white"),
          text = list(size = 6)
        )
      ),
      themes = upset_default_themes(text = element_text(size = 25))
    )
  })
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
