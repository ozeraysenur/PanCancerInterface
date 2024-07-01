library(readr)

# clean_data function : removing unnecessary elements from the data(like html tags)
clean_data <- function(data) {
  
  html_tag_pattern <- "<[^<]+?>"
  
  
  cleaned_data <- gsub(html_tag_pattern, "", data)
  
  cleaned_data <- gsub("\t", " ", cleaned_data)
  
  return(cleaned_data)
}

# load_all_cancer_data function : loading the maftools data to use
load_all_cancer_data <- function() {
  all_data <- list()
  
  cancer_types <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA", 
                    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC", 
                    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ", 
                    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
  
  for (cancer_type in cancer_types) {
    filename <- paste0(cancer_type, "_maftools.maf_modified")
    filepath <- file.path("data", cancer_type, "SNV", filename)
    
    if (file.exists(filepath)) {
      
      data <- readLines(filepath)
      
      
      cleaned_data <- clean_data(data)
      
      
      all_data[[cancer_type]] <- cleaned_data
    } else {
      warning(paste("File not found:", filepath))
    }
  }
  
  return(all_data)
}

# convert_to_matrix function: converting modified data into a numeric matrix
convert_to_matrix <- function(cleaned_data) {
  
  split_data <- strsplit(cleaned_data, " ")
  
  
  numeric_data <- lapply(split_data, function(x) {
    numeric_values <- suppressWarnings(as.numeric(x[-1]))
    # Checking all the values are numeric
    numeric_values[!is.na(numeric_values)]
  })
  
  # Ensure all rows have the same length 
  max_length <- max(sapply(numeric_data, length))
  padded_data <- lapply(numeric_data, function(x) {
    length(x) <- max_length
    return(x)
  })
  
  
  matrix_data <- do.call(rbind, padded_data)
  
  
  rownames(matrix_data) <- sapply(split_data, function(x) gsub('""', '', x[1]))
  
  
  matrix_data[is.na(matrix_data)] <- 0
  
  return(matrix_data)
}


all_cancer_data <- load_all_cancer_data()


cancer_matrices <- list()

# Converting data to a matrix and store it in the list
for (cancer_type in names(all_cancer_data)) {
  cancer_matrices[[cancer_type]] <- convert_to_matrix(all_cancer_data[[cancer_type]])
}

combined_data <- do.call(rbind, lapply(names(cancer_matrices), function(cancer_type) {
  data.frame(Gene = rownames(cancer_matrices[[cancer_type]]),
             CancerType = cancer_type,
             Mutations = rowSums(cancer_matrices[[cancer_type]]))
}))


mutation_percentages <- aggregate(Mutations ~ Gene + CancerType, combined_data, sum)
mutation_percentages <- reshape(mutation_percentages, timevar = "CancerType", idvar = "Gene", direction = "wide")


mutation_matrix <- as.matrix(mutation_percentages[, -1])
rownames(mutation_matrix) <- mutation_percentages$Gene

# replacing NA values with 0s
mutation_matrix[is.na(mutation_matrix)] <- 0

# Create gene pool for gene options
sng_genes <- c("ABCB1", "ACE", "ADRB2", "AGT", "APOB", "AR", "ATM", "CYP1A1", "CYP2C9", "CYP2D6", "CYP3A4", "EGFR", "ERCC2", "ESR1", "FTO", "GSTM1", "GSTP1", "GSTT1", "HLA-B", "IL6", "JAK2", "KRAS", "LCT", "MTHFR", "NAT2", "NOS3", "OPRM1", "PTEN", "SLC6A4", "SLCO1B1", "TGFB1", "TNF", "VDR", "XRCC1")
cnv_genes <- c("BRCA1", "BRCA2", "EGFR", "MYC", "PTEN", "RB1", "TP53", "ERBB2", "CCND1", "CDK4", "CCNE1", "FGFR1", "MDM2", "MET", "PDGFRA", "PIK3CA", "AKT1", "AKT2", "CDKN2A", "CDKN2B", "CTNNB1", "FLT3", "KRAS", "NRAS", "PTCH1", "SMO", "SOX2", "TERT", "VEGFA", "YAP1")
transcriptome_genes <- c("GAPDH", "ACTB", "TPM1", "EGFR", "BRCA1", "BRCA2", "ESR1", "HER2", "MTOR", "MYC", "BCL2", "BAX", "CDH1", "EPCAM", "FOXP3", "GATA3", "HIF1A", "JUN", "MAPK1", "MAPK3", "NF1", "NF2", "PIK3R1", "SMAD4", "STAT3", "TP63", "VEGFA", "VIM", "WT1", "ZEB1")

sng_matrix <- mutation_matrix[rownames(mutation_matrix) %in% sng_genes, , drop = FALSE]
cnv_matrix <- mutation_matrix[rownames(mutation_matrix) %in% cnv_genes, , drop = FALSE]
transcriptome_matrix <- mutation_matrix[rownames(mutation_matrix) %in% transcriptome_genes, , drop = FALSE]
