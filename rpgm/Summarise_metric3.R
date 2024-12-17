
setwd("C:/Users/divya_vq9ublx/Spatial-Metrics-Toolkit")

# Summarize each sample for positive marker frequency and export a --------



# Load required libraries
library(data.table)

# List all .csv.gz files
files <- list.files("data/per-cell", pattern = "\\.csv\\.gz$", full.names = TRUE)

# Initialize an empty list to store results
results <- list()


# Loop through each file
for (file in files) {
  # Read the data using read.csv() with gzfile()
  data <- read.csv(gzfile(file))
  
  # Extract sample name from file name
  sample_name <- tools::file_path_sans_ext(basename(file))
  sample_name <- tools::file_path_sans_ext(sample_name) # Remove .gz
  
  # Identify marker columns
  marker_columns <- grep("CD", colnames(data), value = TRUE)
  
  # Calculate the number of positive cells for each marker
  num_positive <- colSums(data[, marker_columns], na.rm = TRUE)
  
  # Calculate positivity (proportion of positive cells)
  positivity <- num_positive / nrow(data)
  
  
  sample_result <- data.frame(
    Sample = sample_name,
    Marker = names(num_positive),
    NumPositiveCells = num_positive,
    Positivity = positivity
  )
  
  
  results[[sample_name]] <- sample_result
}

# Combine all results into a single data frame
summary_df <- do.call(rbind, results)

# Ensure the output directory exists
if (!dir.exists("outputs/metrics")) {
  dir.create("outputs/metrics", recursive = TRUE)
}

# Save the summary to a CSV file
output_path <- "outputs/metrics/marker_summary3.csv"
write.csv(summary_df, output_path, row.names = FALSE)

cat("Summary saved to:", output_path, "\n")
