library(shiny)
library(data.table)
library(doParallel)
library(reader)
library(ComplexUpset)
library(ggplot2)
library(ggrepel)
library(zip)
library(dplyr)
library(VennDiagram)
library(shinyjs)


ui <- fluidPage(
  titlePanel("FireMAGE Data Visualization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload CSV Files or ZIP file containing CSVs",
                multiple = TRUE,
                accept = c('text/csv',
                           'text/comma-separated-values,text/plain',
                           '.csv',
                           '.zip')),
      uiOutput("species_select"), 
      uiOutput("trait_select"),
      width = 2
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Data Visualization",
                 tabsetPanel(
                   tabPanel("Violin Plot", plotOutput("candidatePlot")),
                   tabPanel("Upset Plot", plotOutput("upsetPlot")),
                   tabPanel("Summary", 
                            div(class = "row",
                                div(class = "col-sm-6", h2("Summary")),
                                div(class = "col-sm-6", style = "text-align: right;", downloadButton("downloadSummary", "Download Summary"))
                            ),
                            dataTableOutput("summaryTable")),
                   tabPanel("Full Table", 
                            div(class = "row",
                                div(class = "col-sm-6", h2("Full table")),
                                div(class = "col-sm-6", style = "text-align: right;", downloadButton("downloadData", "Download Full Table"))
                            ),
                            dataTableOutput("ViewList"))
                 )
        ),
        tabPanel("Multi-trait",
                 tabsetPanel(
                   tabPanel("Genes", 
                            div(class = "row",
                                div(class = "col-sm-6", h2("Multi-trait Genes")),
                                div(class = "col-sm-6", style = "text-align: right;", downloadButton("downloadMultiTraitGenes", "Download Multi-trait Genes"))
                            ),
                            uiOutput("trait_select_gene"), 
                            dataTableOutput("uniqueGeneNamesTable")),
                   tabPanel("Orthogroups", 
                            div(class = "row",
                                div(class = "col-sm-6", h2("Multi-trait Orthogroups")),
                                div(class = "col-sm-6", style = "text-align: right;", downloadButton("downloadMultiTraitOrthogroups", "Download Multi-trait Orthogroups"))
                            ),
                            dataTableOutput("multiTraitOrthogroupsTable")),
                   tabPanel("Venn Diagram Visualization",
                            div(class = "row",
                                div(class = "col-sm-6", h2("Venn Diagram Visualization")),
                                div(class = "col-sm-6", style = "text-align: right;", downloadButton("downloadVennDiagram", "Download Multi-trait Venn Diagram"))
                            ),
                            plotOutput("vennPlot")
                            )
                 )
        ),
        tabPanel("Help",
                 fluidRow(
                   column(width = 8, 
                          div(class = "row", 
                              h3("General"),
                              p("This is a tool for analyzing genetic data related to various species and traits."),
                              h3("Input"),
                              p("Upload CSV files or a ZIP file containing CSVs. Each CSV file should contain the following columns: org, loci, Orthogroup, trait, Gene Name, Description, present, genecount, pFDR."),
                              h3("Species Selection"),
                              p("Select the species to be included in the analysis from the dropdown list. Multiple species can be selected."),
                              h3("Trait Selection"),
                              p("Select the traits to be included in the analysis from the dropdown list. Multiple traits can be selected."),
                              h3("Plots"),
                              p("Violin Plot: This plot visualizes the distribution of FDR (False Discovery Rate) values for different traits and levels of species presence. The width of the 'violin' shape at any given y value (FDR) represents the proportion of data that falls at that level. It's a useful plot for understanding the distribution of your data and spotting any outliers or significant clusters."),
                              p("Upset Plot: This plot displays the intersection of traits for different species. It is similar to a Venn diagram but more suitable for datasets with many sets and complex overlaps. Each row of the matrix at the bottom represents a set (in this case, a species), and each column of the matrix represents a combination of intersections. The bar chart above the matrix indicates the size of each intersection."),
                              h3("Tables"),
                              p("Summary Table: View a summary of the gene count, mean FDR, and standard deviation of FDR for each trait and species."),
                              p("Full Table: View the complete data for the selected species and traits."),
                              p("Multi-trait Genes: View the genes associated with multiple traits."),
                              h3("Download"),
                              p("You can download the Summary Table, Full Table, and Multi-trait Genes Table as CSV files.")
                          )
                   )
                 )
        )
      )
    )
  )
)

# Define server logic 
server <- function(input, output, session) {
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
  
  fullTable <- reactiveVal()  # Create a reactive variable to hold the full table
  summaryTable <- reactiveVal()  # Create a reactive variable to hold the summary table
  
  output$ViewList <- renderDataTable({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
    fullTable(table)  # Update the reactive variable for the full table
    table
  })
  
  output$summaryTable <- renderDataTable({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
    summary_stats <- table[, .(Genes = length(unique(`Gene Name`)), Mean_genes_per_loci = mean(genecount), SD_genes_per_loci = sd(genecount)), by = .(trait,org,present)]
    summaryTable(summary_stats)  # Update the reactive variable for the summary table
    summary_stats
  })
  
  # Download full table
  output$downloadData <- downloadHandler(
    filename = "full_table.csv",
    content = function(file) {
      write.csv(fullTable(), file, row.names = FALSE)
    }
  )
  
  # Download summary table
  output$downloadSummary <- downloadHandler(
    filename = "summary_table.csv",
    content = function(file) {
      write.csv(summaryTable(), file, row.names = FALSE)
    }
  )
  
  # Creating the dropdown for trait selection in Gene Names
  output$trait_select_gene <- renderUI({
    if (is.null(candidateList())) return(NULL)
    candidateList <- candidateList()
    selectInput("selected_trait_gene", "Select Trait for Gene Names", choices = unique(candidateList$trait), selected = unique(candidateList$trait)[1], multiple = TRUE)
  })
  
  
  # Adjust this data table for unique gene names
  uniqueGeneNamesTable <- reactiveVal()  # Create a reactive variable to hold the unique gene names table
  
  output$uniqueGeneNamesTable <- renderDataTable({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait_gene)) return()
    table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait_gene]
    unique_gene_table <- table[, .(Traits = toString(unique(trait))), by = `Gene Name`]
    uniqueGeneNamesTable(unique_gene_table)  # Update the reactive variable for the unique gene names table
    unique_gene_table
  })
  
  output$uniqueGeneNamesTable <- renderDataTable({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait_gene)) return()
    table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait_gene]
    
    # Group by 'Gene Name', 'Orthogroup' and 'Description', then aggregate the traits
    unique_gene_table <- table[, .(Traits = toString(unique(trait)), Orthogroup = first(Orthogroup), Description = first(Description)), by = `Gene Name`]
    
    # Filter for only multi-trait genes and add count of traits
    unique_gene_table <- unique_gene_table[unlist(lapply(strsplit(unique_gene_table$Traits, split = ", "), length)) > 1]
    
    # Add a new column for the count of traits
    unique_gene_table[, 'TraitCount' := unlist(lapply(strsplit(Traits, split = ", "), length))]
    
    # Order by TraitCount in descending order
    unique_gene_table <- unique_gene_table[order(-TraitCount)]
    
    uniqueGeneNamesTable(unique_gene_table)  # Update the reactive variable for the unique gene names table
    unique_gene_table
  })
  
  
  
  
  # Download multi-trait genes table
  output$downloadMultiTraitGenes <- downloadHandler(
    filename = "multi_trait_genes_table.csv",
    content = function(file) {
      write.csv(uniqueGeneNamesTable(), file, row.names = FALSE)
    }
  )
  
  # Create a reactive variable to hold the multi-trait orthogroups table
  multiTraitOrthogroupsTable <- reactiveVal() 
  
  # In the server function
  output$multiTraitOrthogroupsTable <- renderDataTable({
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
    
    # Group by 'Orthogroup', and aggregate the traits
    multi_trait_orthogroups <- table[, .(Traits = toString(unique(trait))), by = Orthogroup]
    
    # Filter for only multi-trait orthogroups and add count of traits
    multi_trait_orthogroups <- multi_trait_orthogroups[unlist(lapply(strsplit(multi_trait_orthogroups$Traits, split = ", "), length)) > 1]
    
    # Add a new column for the count of traits
    multi_trait_orthogroups[, 'TraitCount' := unlist(lapply(strsplit(Traits, split = ", "), length))]
    
    # Now, create the species present/absent columns
    species_list <- unique(table$org)  # Get the unique species from the data
    for (species in species_list) {
      multi_trait_orthogroups[, paste0("IsPresent_", species) := Orthogroup %in% table[org == species]$Orthogroup]
    }
    
    # Order by TraitCount in descending order
    multi_trait_orthogroups <- multi_trait_orthogroups[order(-TraitCount)]
    
    multiTraitOrthogroupsTable(multi_trait_orthogroups)  # Update the reactive variable for the multi-trait orthogroups table
    multi_trait_orthogroups
  })
  
  # In the download button function
  output$downloadMultiTraitOrthogroups <- downloadHandler(
    filename = "multi_trait_orthogroups_table.csv",
    content = function(file) {
      write.csv(multiTraitOrthogroupsTable(), file, row.names = FALSE)
    }
  )
  
  display_venn <- function(x, ...){
    library(VennDiagram)
    grid.newpage()
    venn_object <- venn.diagram(x, filename = NULL, ...)
    grid.draw(venn_object)
  }
  
  output$vennPlot <- renderPlot({
    
    if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
    
    selected_orthogroup <- "OG0000141"
    
    filtered_data <- candidateList()[trait %in% input$selected_trait & Orthogroup == selected_orthogroup]
    if (nrow(filtered_data) == 0) return(NULL)
    
    venn_traits <- selected_traits()
    
    venn_data <- lapply(venn_traits, function(trait_name){
      filtered_data[trait == trait_name, unique(`Gene Name`)]
    })
    
    names(venn_data) <- venn_traits
    
    display_venn(venn_data,
      category.names = names(venn_data)
    )
  })
  # output$multiTraitVennDiagram <- renderDataTable({
  #   if (is.null(candidateList()) | is.null(input$selected_species) | is.null(input$selected_trait)) return()
  #   table <- candidateList()[org %in% input$selected_species & trait %in% input$selected_trait]
  # }
  # )
}

# Run the application 
shinyApp(ui = ui, server = server)
