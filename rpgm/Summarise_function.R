# Load required libraries
library(data.table)
# 
# # Function to summarize marker positivity
# generate_marker_summary <- function(input_folder, output_folder) {
#   # List all .csv.gz files
#   files <- list.files(input_folder, pattern = "\\.csv\\.gz$", full.names = TRUE)
#   
#   # Initialize an empty list to store results
#   results <- list()
#   
#   # Loop through each file
#   for (file in files) {
#     # Read the data using read.csv() with gzfile()
#     data <- read.csv(gzfile(file))
#     
#     # Extract sample name from file name
#     sample_name <- tools::file_path_sans_ext(basename(file))
#     sample_name <- tools::file_path_sans_ext(sample_name) # Remove .gz
#     
#     # Identify marker columns
#     marker_columns <- grep("CD", colnames(data), value = TRUE)
#     
#     # Calculate the number of positive cells for each marker
#     num_positive <- colSums(data[, marker_columns], na.rm = TRUE)
#     
#     # Calculate positivity (proportion of positive cells)
#     positivity <- num_positive / nrow(data)
#     
#     # Combine results into a data frame
#     sample_result <- data.frame(
#       Sample = sample_name,
#       Marker = names(num_positive),
#       NumPositiveCells = num_positive,
#       Positivity = positivity
#     )
#     
#     # Add the sample data frame to the results list
#     results[[sample_name]] <- sample_result
#   }
#   
#   # Combine all results into a single data frame
#   summary_df <- do.call(rbind, results)
#   
#   # Ensure the output directory exists
#   if (!dir.exists(output_folder)) {
#     dir.create(output_folder, recursive = TRUE)
#   }
#   
#   # Save the summary to a CSV file
#   output_path <- file.path(output_folder, "marker_summary.csv")
#   write.csv(summary_df, output_path, row.names = FALSE)
#   
#   cat("Summary saved to:", output_path, "\n")
# }
# 
# # Example usage
# generate_marker_summary("data/per-cell", "outputs/metrics")
library(dplyr)#summarise functionality
#try a cleaner way perhaps?
generate_marker_summary = function(config, path){
  df = fread(path, data.table = FALSE)
  #if there is a tissue classifier
  if(!is.null(config$variables$tissue_class_label)){
    overall = summarize_markers(config, df) %>%
      mutate(!!config$variables$tissue_class_label := "Overall")
    df_list = split(df, df[[config$variables$tissue_class_label]])
    df_summs = lapply(df_list, function(df1){
      summarize_markers(config, df1)
    }) %>%
      do.call(bind_rows, .) %>%
      mutate(!!config$variables$tissue_class_label := names(df_list)) %>%
      bind_rows(overall, .)
  } else {
    df_summs = summarize_markers(config, df)
  }
  
  return(df_summs)
}

summarize_markers = function(config, df){
  df %>%
    group_by(get(config$variables$sample_id)) %>%
    summarise(across(c(!!config$variables$markers),
                     ~ sum(.x),
                     .names = "{col} Cells"),
              `Total Cells` = n()) %>%
    mutate(across(contains(!!config$variables$markers),
                  ~ .x / `Total Cells` * 100,
                  .names = "% {col}"),
           .before  = `Total Cells`)
}
