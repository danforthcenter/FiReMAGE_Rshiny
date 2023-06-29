library(shiny)
library(readr)
library(ggplot2)
library(plotly)
library(dplyr)
library(RColorBrewer)
library(UpSetR)
library(tidyr)

# User interface
# User interface
ui <- fluidPage(
  tags$head(
    tags$style(
      type = "text/css",
      ".scrollable_plot {
        height: 700px;
        overflow-y: scroll;
        overflow-x: hidden;
      }"
    )
  ),
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Please choose a CSV file",
                multiple = TRUE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      uiOutput("orgInput"), # Dynamic UI for organisms
      uiOutput("traitInput"), # Dynamic UI for traits
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Scatter Plot", div(class = "scrollable_plot", plotlyOutput("scatterPlot"))),
        tabPanel("Histogram", plotlyOutput("histogramPlot")),
        tabPanel("UpSet Plot", plotOutput("upsetPlot"))
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
  
  # Making big Palette for traits
  bigPalette <- colorRampPalette(brewer.pal(12, "Paired"))
  subPalette <- bigPalette(21)
  
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
  
  
  # Subset data based on user selection
  selectedData <- reactive({
    req(input$org, input$trait)
    data() %>%
      filter(org %in% input$org,
             trait %in% input$trait)
  })


  #Scatter Plot
  output$scatterPlot <- renderPlotly({
    req(selectedData())
    subPalette <- rainbow(length(unique(selectedData()$trait)))
    
    p <- selectedData() %>%
      ggplot(aes(y = genecount, x = trait, color = trait)) +
      geom_boxplot() +
      scale_color_manual(values = subPalette) +
      facet_wrap(~org, scales = "free", ncol = 1, strip.position = "left") +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20), 
        strip.text = element_text(size = 14, face = "bold", angle = 0), # Rotate facet labels to be horizontal
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10, angle = 45, hjust = 1), # Rotate x-axis labels
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)
      ) +
      coord_cartesian(clip = "off") + # prevent clipping of strip labels
      theme(plot.margin = margin(1, 1, 1, 1, "cm")) # Adjust plot margins
    
    # Increase overall plot height; adjust based on the number of facets
    plot_height <- 300 * length(unique(selectedData()$org))
    
    ggplotly(p, height = plot_height, width = 1200) # Increase plot width
  })
  
  output$histogramPlot <- renderPlotly({
    req(selectedData())
    plot_ly(selectedData(), x = ~genecount, type = 'histogram', color = ~trait) %>%
      layout(yaxis = list(title = "Count"),
             xaxis = list(title = "Gene Count"),
             autosize = TRUE)
  })

  # UpSet Plot
  output$upsetPlot <- renderPlot({
    req(selectedData())
    
    # Prepare data for UpSet plot: make a list of sets
    upset_data <- lapply(unique(selectedData()$org), function(org) {
      filter(selectedData(), org == org)$loci
    })
    
    names(upset_data) <- unique(selectedData()$org)
    
    # Make UpSetR plot
    upset(fromList(upset_data), 
          nsets = length(unique(selectedData()$org)), 
          nintersects = 30, 
          order.by = "freq", 
          mainbar.y.label = "Gene count", 
          sets.x.label = "Organism count", 
          text.scale = c(1.3, 1.3, 1, 1, 2, 0.7))
  })
  
}
# Run the application 
shinyApp(ui = ui, server = server)
