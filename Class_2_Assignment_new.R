#HOPE EDIGHE JULIUS_Class_2
#FOR INSTRUCTOR CONSIDERATION-summary_table was mistakenly omitted in the RData uploaded in the google form  already submitted
#ASSIGNMENT
# Define the input folder (where raw data files are stored) and the output folder (where results will be saved).

input_dir <- "DEGs_Raw_Data"

output_dir <- "DEGs_Result"

# create output folder if not already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# create input folder if not already exist

if (!dir.exists(input_dir)) {
  dir.create(input_dir)
}

# Write a function classify_gene 
classify_genes <- function(logFC, padj) {
  ifelse(logFC > 1 & padj < 0.05, "Upregulated",
         ifelse(logFC < -1 & padj < 0.05, "Downregulated",
                "Not_Significant"))
}

# List which files to process
file_to_process <- c("DEGs_data_1.csv", "DEGs_data_2.csv")

result_list <- list()

#Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)

for (file_names in file_to_process) {
  cat("\nProcessing:", file_names, "\n")
  
  input_file_path <- file.path(input_dir, file_names)
  
#Import dataset  
  data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. checking for missing values...\n") 
  
#Replace missing padj values with 1
  if("padj" %in% names(data)){
    missing_count <- sum(is.na(data$padj))
    
    cat("Missing values in 'padj':", missing_count, "\n")
    data$padj[is.na(data$padj)] <- 1
  }
  
    data$classification <- classify_genes(data$logFC, data$padj)
    cat(" Genes has been classified successfully.\n")
    
#Add a new column 'status'
    data$status <- classify_genes(data$logFC, data$padj)
    cat("Status column has been created.\n")

#Save processed files into Results folder    
    result_list[[file_names]] <- data
  
    output_file_path <- file.path(output_dir, paste0("DEGs_result", file_names))
    write.csv(data, output_file_path, row.names =  FALSE)
    cat("Results saved to:", output_file_path, "\n")
}

results_1 <- result_list[[1]] 
results_2 <- result_list[[2]]

for (file_name in file_to_process) {
  cat("\nSummary for:", file_name, "\n")
  
#Retrieve the data frame for the current file from the result_list
  data <- result_list[[file_name]]
  
#Generate and print the summary table for the 'classification' column
  summary_table <- table(data$classification)
  print(summary_table)
}

#FOR INSTRUCTOR CONSIDERATION-summary_table was mistakenly omitted in the RData uploaded in the google form  already submitted