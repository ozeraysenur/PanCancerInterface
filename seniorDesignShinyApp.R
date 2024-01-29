

library(shiny)


ui <- fluidPage(

  
    titlePanel("TCGAnalyzeR Pan-Cancer View"),

    sidebarLayout(
        selectInput("myGenes",
                    "My Genes",
                    choices = c("Single Nucleotid", "Copy Number Var", 
                                "Transcriptome"),
                    selected = "Single Nucleotid"
                    
          
        ),

        
        mainPanel(
           plotlyOutput("genePlot")
        )
    )
)

library(plotly)

server <- function(input, output) {
  
  output$genePlot <- renderPlotly({
    
    genes <- combined_paad_splitted$external_gene_name
    
    if (input$myGenes == "Single Nucleotid") {
      y <- combined_paad_splitted$t
    } else if (input$myGenes == "Copy Number Var") {
      y <- combined_paad_splitted$P.Value
    } else if (input$myGenes == "Transcriptome") {
      y <- combined_paad_splitted$adj.P.Val
    }
    
    plot_ly(x = ~genes, y = ~y, type = 'scatter', mode = 'markers', marker = list(color = 'pink')) 
    
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
