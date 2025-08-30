#HOPE EDIGHE JULIUS_Class_2

#ASSIGNMENT

input_dir <- "DEGs_Raw_Data"

output_dir <- "DEGs_Result"


if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

if (!dir.exists(input_dir)) {
  dir.create(input_dir)
}

classify_genes <- function(logFC, padj) {
  ifelse(logFC > 1 & padj < 0.05, "Upregulated",
         ifelse(logFC < -1 & padj < 0.05, "Downregulated",
                "Not_Significant"))
}


file_to_process <- c("DEGs_data_1.csv", "DEGs_data_2.csv")

result_list <- list()


for (file_names in file_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. checking for missing values...\n") 
  
  
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
    data$classification <- classify_genes(data$logFC, data$padj)
    cat(" Genes has been classified successfully.\n")
    
    data$status <- classify_genes(data$logFC, data$padj)
    cat("Status column has been created.\n")
    
    result_list[[file_names]] <- data
  
    output_file_path <- file.path(output_dir, paste0("DEGs_result", file_names))
    write.csv(data, output_file_path, row.names =  FALSE)
    cat("Results saved to:", output_file_path, "\n")
}

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]
