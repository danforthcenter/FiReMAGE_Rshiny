library(shiny)
library(readr)
library(plotly)
library(shinydashboard)
library(reactable)
library(shinyWidgets)
library(RColorBrewer)

# Load the data outside of the server function
snps_subset <- read_csv("C:/Users/srava/Desktop/Intern/snps_subset.csv")

# Get unique species from the data
unique_species <- unique(snps_subset$org)

# Define UI
ui <- dashboardPage(
  dashboardHeader(title = "SNPs Subset Viewer"),
  dashboardSidebar(
    selectInput("sort", "Sort by:", choices = colnames(snps_subset)),
    numericInput("rows", "Number of rows to display:", min = 1, max = nrow(snps_subset), value = 10),
    pickerInput("filter_trait", "Filter by Trait:", choices = unique(snps_subset$trait), multiple = TRUE, options = list(`actions-box` = TRUE)),
    pickerInput("filter_species", "Filter by Species:", choices = unique_species, multiple = TRUE, options = list(`actions-box` = TRUE)),
    selectInput("chromosome", "Filter by Chromosome:", choices = unique(snps_subset$chr))
  ),
  dashboardBody(
    tabsetPanel(
      tabPanel("Table", reactableOutput("table")),
      tabPanel("Scatter Plots", uiOutput("scatterPlotTabs")),
      tabPanel("Histogram", plotlyOutput("histogram"))
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Reactive expression to create a subset and sort
  dataSubset <- reactive({
    data <- snps_subset
    
    # Filter by trait if requested
    if (!is.null(input$filter_trait) && length(input$filter_trait) > 0) {
      data <- data[data$trait %in% input$filter_trait,]
    }
    
    # Filter by species if requested
    if (!is.null(input$filter_species) && length(input$filter_species) > 0) {
      data <- data[data$org %in% input$filter_species,]
    }
    
    # Filter by chromosome if requested
    if (input$chromosome != "All") {
      data <- data[data$chr == input$chromosome,]
    }
    
    # Sort the data
    data <- data[order(data[[input$sort]]),]
    
    # Limit the number of rows
    data <- head(data, input$rows)
    
    data
  })
  
  # Create the table output object
  output$table <- renderReactable({
    reactable(dataSubset())
  })
  
  # Create individual scatter plots for each species
  output$scatterPlotTabs <- renderUI({
    lapply(unique_species, function(species) {
      plotlyOutput(paste0("scatterPlot_", species))
    })
  })
  
  observe({
    lapply(unique_species, function(species) {
      output[[paste0("scatterPlot_", species)]] <- renderPlotly({
        filtered_data <- dataSubset()
        filtered_data <- filtered_data[filtered_data$org == species, ]
        
        plot_ly(filtered_data, x = ~chr, y = ~bp, type = 'scatter', mode = 'markers', text = ~trait, hoverinfo = 'text') %>%
          layout(xaxis = list(title = "Chromosome"), yaxis = list(title = "BP"))
      })
    })
  })
  
  # Create histogram
  output$histogram <- renderPlotly({
    plot_ly(dataSubset(), x = ~genecount, type = 'histogram')
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
