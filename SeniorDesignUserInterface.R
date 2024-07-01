library(shiny)
library(pheatmap)

sng_genes <- c("ABCB1", "ACE", "ADRB2", "AGT", "APOB", "AR", "ATM", "CYP1A1", "CYP2C9", "CYP2D6", "CYP3A4", "EGFR", "ERCC2", "ESR1", "FTO", "GSTM1", "GSTP1", "GSTT1", "HLA-B", "IL6", "JAK2", "KRAS", "LCT", "MTHFR", "NAT2", "NOS3", "OPRM1", "PTEN", "SLC6A4", "SLCO1B1", "TGFB1", "TNF", "VDR", "XRCC1", "TP53")
cnv_genes <- c("BRCA1", "BRCA2", "EGFR", "MYC", "PTEN", "RB1", "TP53", "ERBB2", "CCND1", "CDK4", "CCNE1", "FGFR1", "MDM2", "MET", "PDGFRA", "PIK3CA", "AKT1", "AKT2", "CDKN2A", "CDKN2B", "CTNNB1", "FLT3", "KRAS", "NRAS", "PTCH1", "SMO", "SOX2", "TERT", "VEGFA", "YAP1", "TP53")
transcriptome_genes <- c("GAPDH", "ACTB", "TPM1", "EGFR", "BRCA1", "BRCA2", "ESR1", "HER2", "MTOR", "MYC", "BCL2", "BAX", "CDH1", "EPCAM", "FOXP3", "GATA3", "HIF1A", "JUN", "MAPK1", "MAPK3", "NF1", "NF2", "PIK3R1", "SMAD4", "STAT3", "TP63", "VEGFA", "VIM", "WT1", "ZEB1", "TP53")
C:/Users/ozera/OneDrive/Masa??st??/SeniorDesignUserInterface.R

# Remove the "Mutations" prefix from the column names
colnames(mutation_matrix) <- gsub("Mutations\\.", "", colnames(mutation_matrix))

sng_matrix <- mutation_matrix[rownames(mutation_matrix) %in% sng_genes, , drop = FALSE]
cnv_matrix <- mutation_matrix[rownames(mutation_matrix) %in% cnv_genes, , drop = FALSE]
transcriptome_matrix <- mutation_matrix[rownames(mutation_matrix) %in% transcriptome_genes, , drop = FALSE]


set.seed(123) 
cancer_types <- colnames(mutation_matrix)
msi_high <- sample(0:5, length(cancer_types), replace = TRUE)
names(msi_high) <- cancer_types
msi_low <- sample(0:15, length(cancer_types), replace = TRUE)
names(msi_low) <- cancer_types
immune_positive <- sample(0:10, length(cancer_types), replace = TRUE)
names(immune_positive) <- cancer_types
immune_negative <- sample(0:20, length(cancer_types), replace = TRUE)
names(immune_negative) <- cancer_types


get_filtered_matrix <- function(filter_type, filter_value) {
  if (filter_type == "MSI") {
    if (filter_value == "high") {
      msi_matrix <- matrix(msi_high, nrow = 1)
      rownames(msi_matrix) <- "MSI-High"
    } else {
      msi_matrix <- matrix(msi_low, nrow = 1)
      rownames(msi_matrix) <- "MSI-Low"
    }
    colnames(msi_matrix) <- cancer_types
    return(msi_matrix)
  } else if (filter_type == "Immune") {
    if (filter_value == "positive") {
      immune_matrix <- matrix(immune_positive, nrow = 1)
      rownames(immune_matrix) <- "Immune-Positive"
    } else {
      immune_matrix <- matrix(immune_negative, nrow = 1)
      rownames(immune_matrix) <- "Immune-Negative"
    }
    colnames(immune_matrix) <- cancer_types
    return(immune_matrix)
  } else {
    return(NULL)
  }
}

get_heatmap_colors <- function(heatmap_type) {
  if (heatmap_type == "sng") {
    return(colorRampPalette(c("beige", "orange", "#023E8A"))(100))
  } else if (heatmap_type == "cnv") {
    return(colorRampPalette(c("beige", "#8DA750", "purple"))(100))
  } else if (heatmap_type == "transcriptome") {
    return(colorRampPalette(c("#FFF394", "white", "brown"))(100))
  } else {
    return(colorRampPalette(c("gray", "white", "black"))(100))
  }
}

shinyApp(
  ui = fluidPage(
    tags$head(
      tags$style(HTML("
        .heatmap-title {
          font-size: 24px;
          font-weight: bold;
          margin-bottom: 20px;
          text-align: center;
        }
      "))
    ),
    titlePanel("TCGAnalyzeR Pan-Cancer View"),
    sidebarLayout(
      sidebarPanel(
        selectInput("gene_type", "Gene Type",
                    choices = c("Single Nucleotide" = "sng",
                                "Copy Number Var" = "cnv",
                                "Transcriptome" = "transcriptome")),
        selectizeInput("genes", "MyGenes", choices = NULL, multiple = TRUE),
        actionButton("add_gene", "+ Add"),
        selectInput("filter_by", "Filter By",
                    choices = c("None", "MSI", "Immune")),
        uiOutput("filter_type_ui")
      ),
      mainPanel(
        h3(textOutput("heatmap_title"), class = "heatmap-title"),
        plotOutput("heatmap_plot", height = "600px")
      )
    )
  ),
  
  server = function(input, output, session) {
    
    observe({
      updateSelectizeInput(session, "genes", choices = rownames(mutation_matrix), server = TRUE)
    })
    
    observeEvent(input$filter_by, {
      output$filter_type_ui <- renderUI({
        if (input$filter_by == "MSI") {
          selectInput("filter_type", "MSI Type",
                      choices = c("MSI-High" = "high", "MSI-Low" = "low"))
        } else if (input$filter_by == "Immune") {
          selectInput("filter_type", "Immune Type",
                      choices = c("Immune-Positive" = "positive", "Immune-Negative" = "negative"))
        } else {
          return(NULL)
        }
      })
    })
    
    generate_heatmap <- function(heatmap_type, selected_genes, filter_by, filter_type) {
      if (filter_by %in% c("MSI", "Immune")) {
        matrix_data <- get_filtered_matrix(filter_by, filter_type)
      } else {
        matrix_data <- switch(heatmap_type,
                              "sng" = sng_matrix,
                              "cnv" = cnv_matrix,
                              "transcriptome" = transcriptome_matrix)
        
        if (length(selected_genes) > 0) {
          matrix_data <- matrix_data[rownames(matrix_data) %in% selected_genes, , drop = FALSE]
        }
      }
      
      colors <- get_heatmap_colors(heatmap_type)
      pheatmap(matrix_data, cluster_rows = FALSE, cluster_cols = FALSE, 
               show_rownames = TRUE, show_colnames = TRUE, color = colors)
    }
    
    output$heatmap_title <- renderText({
      if (input$filter_by == "MSI") {
        paste("Mutation Number of", ifelse(input$filter_type == "high", "MSI-High", "MSI-Low"), "Gene Regulation Heatmap")
      } else if (input$filter_by == "Immune") {
        paste("Mutation Number of", ifelse(input$filter_type == "positive", "Immune-Positive", "Immune-Negative"), "Gene Regulation Heatmap")
      } else {
        switch(input$gene_type,
               "sng" = "Mutation Number of Single Nucleotide Gene Regulation Heatmap",
               "cnv" = "Mutation Number of Copy Number Variations Gene Regulation Heatmap",
               "transcriptome" = "Mutation Number of Transcriptome Gene Regulation Heatmap")
      }
    })
    
    output$heatmap_plot <- renderPlot({
      generate_heatmap(input$gene_type, input$genes, input$filter_by, input$filter_type)
    })
  }
)
