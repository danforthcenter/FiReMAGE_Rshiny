library(shiny)
library(readr)
library(ggplot2)
library(plotly)
library(dplyr)

# User interface
ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Please choose a CSV file",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      uiOutput("orgInput"), # Dynamic UI for organisms
      uiOutput("traitInput"), # Dynamic UI for traits
      uiOutput("chrInput") # Dynamic UI for chromosomes
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot", plotlyOutput("scatterPlot")),
        tabPanel("Histogram", plotlyOutput("histogramPlot"))
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  
  # Reactive expression to read in the data when uploaded
  data <- reactive({
    req(input$file1)
    read_csv(input$file1$datapath)
  })
  
  # Dynamic UI for organism selection
  output$orgInput <- renderUI({
    selectInput("org", "Select Organism(s)",
                choices = unique(data()$org), multiple = TRUE)
  })
  
  # Dynamic UI for trait selection
  output$traitInput <- renderUI({
    selectInput("trait", "Select Trait(s)",
                choices = unique(data()$trait), multiple = TRUE)
  })
  
  # Dynamic UI for chromosome selection
  output$chrInput <- renderUI({
    selectInput("chr", "Select Chromosome(s)",
                choices = unique(data()$chr), multiple = TRUE)
  })
  
  # Subset data based on user selection
  selectedData <- reactive({
    req(input$org, input$trait, input$chr)
    data() %>%
      filter(org %in% input$org,
             trait %in% input$trait,
             chr %in% input$chr)
  })
  
  # Scatter Plot
  output$scatterPlot <- renderPlotly({
    req(selectedData())
    plot_ly(data = selectedData(), 
            x = ~bp, y = ~genecount, 
            type = 'scatter', mode = 'markers', 
            color = ~trait)
  })
  
  # Histogram
  output$histogramPlot <- renderPlotly({
    req(selectedData())
    plot_ly(data = selectedData(), 
            x = ~bp, 
            type = 'histogram', 
            color = ~trait)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

