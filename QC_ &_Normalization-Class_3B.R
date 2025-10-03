if(!requireNamespace("BiocManager", quietly = TRUE ))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","limma","arrayQualityMetrics",
                       "AnnotationDbi","hgu133plus2.db")) 

install.packages("dplyr")


library(GEOquery)
library(affy)
library(limma)
library(arrayQualityMetrics)
library(AnnotationDbi)
library(hgu133plus2.db)

install.packages("plyr")
library(dplyr)


gse_data <- getGEO(filename = "C:/Users/USER/Desktop/GSE31547/GSE31547_series_matrix.txt.gz")


# Load the local file WITHOUT trying to download the platform file
gse_data <- getGEO(filename = "C:/Users/USER/Desktop/GSE31547/GSE31547_series_matrix.txt.gz", getGPL = FALSE)


# Extract the expression data matrix directly from the ExpressionSet object
expression_data <- exprs(gse_data)
feature_data <- fData(gse_data)

BiocManager::install(c("GEOquery","affy","limma","arrayQualityMetrics","AnnotationDbi","hgu133plus2.db"))


phenotype_data <- pData(gse_data)

sum(is.na(phenotype_data$soure_name_ch1))

# Define the names of the folders you want to create
project_folders <- c("R_Scripts", "Raw_Data", "processed_Data", "Results")

# Loop through each folder name
for (folder in project_folders) {
  # Check if the folder already exists
  if (!dir.exists(folder)) {
    # If it doesn't exist, create it
    dir.create(folder)
    cat("Created folder:", folder, "\n")
  } else {
    # If it exists, print a message
    cat("Folder already exists:", folder, "\n")
  }
}

# Define the source and destination paths
source_path <- "C:/Users/USER/Desktop/GSE31547/GSE31547_RAW.tar"
destination_path <- "Raw_Data/GSE31547_RAW.tar"

source_path <- "C:/Users/USER/Desktop/GSE31547/GSE31547_RAW.tar"
destination_dir <- "Raw_Data"
destination_path <- file.path(destination_dir, basename(source_path))

if (file.exists(source_path)) {
  if (!dir.exists(destination_dir)) {
    dir.create(destination_dir, recursive = TRUE)
  }
  file.rename(from = source_path, to = destination_path)
  cat("Moved", basename(source_path), "to the Raw_Data folder.\n")
} else {
  stop("Source file not found at:", source_path)
}

untar("Raw_Data/GSE31547_RAW.tar", exdir = "Raw_Data/CEL_Files" )
library(affy)

raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files" )
raw_data

#1. QC Before Pre-processing

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_Raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)
##There are 5 outliers (sample array data 6, 18, 28, 40, 50)


#RMA Normalization
normalized_data <- rma(raw_data)

#QC after normalization
arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_Normalized_Data",
                    force = TRUE)
##There are 4 outliers (sample array data 10, 18, 20, 30)


#Removing outliers
sample_id_column <- "sample"

# The names of the outlier samples to remove
outlier_names <- c("10", "18", "20", "30")

# Create a logical vector: TRUE for samples to keep, FALSE for samples to remove
samples_to_keep <- !(pData(normalized_data)[[sample_id_column]] %in% outlier_names)

# Subset the ExpressionSet object to create a clean version
normalized_data_clean <- normalized_data[, samples_to_keep]

# Optional: Verify the removal
dim(normalized_data)
dim(normalized_data_clean)



#2. Apply filtering to remove low-intensity probes

processed_data <- as.data.frame(exprs(normalized_data_clean))
dim(processed_data)

row_median <- rowMedians(as.matrix(processed_data))
row_median

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution")

threshold <- 4


abline(v = threshold, col = "black", lwd = 2)
indx <- row_median > threshold
filtered_data <- processed_data[indx, ]

processed_data <- filtered_data
##Transcripts remaining are 21703

#3.Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)
class(phenotype_data$source_name_ch1)

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("Normal lung", "Lung adenocarcinoma"),
                 label = c("normal", "cancer"))
class(groups)
levels(groups)
