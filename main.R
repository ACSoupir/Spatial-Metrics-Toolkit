#run with
#"Rscript main.R -y config.yml"
#moffitt hpc "module load gcc/11.2.0 " dec 18 2024

#global libraries
cat("Prepping environment\n")
# List of packages
packages <- c("yaml", "optparse", "parallel", "data.table",
              "ggplot2", "tibble", "dplyr", "tidyr", "Polychrome", "ggpubr",
              "spatstat.geom", "spatstat.explore", "dbscan",
              "R.utils", "rslurm")

# Install missing packages 
missing_packages <- setdiff(packages,names(installed.packages()[,1]))
install.packages(missing_packages, quiet = TRUE, repos = "https://cloud.r-project.org")

# Suppress installation messages from printing to console
loaded = lapply(packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
})

#for testing
#option_list = read_yaml("config.yml")
cat("Making options list\n")
#option - only taking in a yaml
option_list = list(
  make_option(
    opt_str = c("-y", "--yaml"),
    type = 'character',
    default = "config.yml",
    help = "path to yaml file with parameters",
    metavar = 'character'
  ),
  make_option(
    opt_str = c("-c", "--cores"),
    type = 'integer',
    default = 1,
    help = "number of cores for parallel processing; config overrides flag",
    metavar = 'character'
  )
)

#parse the option to get the yaml location
cat("Parsing Options\n")
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#read in the yaml
cat("Reading yaml\n")
config = read_yaml(opt$yaml)

# Check if SLURM command 'sinfo' is available
if (Sys.which("sinfo") != "" & config$slurm == TRUE) {
  cat("SLURM is available. 'sinfo' command found.\n")
  config$slurm = TRUE
} else {
  cat("SLURM not found. 'sinfo' command is not available.\n")
  config$slurm = FALSE
}
#source functions
cat("Sourcing functions\n")
source("rpgm/source_functions.R")

#source validation functions
cat("Running validations\n")
validate_config(config)
# example call: function to validate df columns
# validate_df_columns(df, c("col1", "col2"))

#if missing cores in config replace
if(is.null(config$cores)){
  config$cores = opt$cores
}

#spatial files
cat("Identifying spatial files\n")
sfiles = list.files(config$paths$spatial, 
                    pattern = 'csv',
                    full.names = TRUE)

#prep folder structure
cat("Creating folder structure\n")
createfolders(config)

#plot functions for raw csv inputs
#use data.table to import gz
cat("Plotting xy points\n")
if(config$slurm){
  tmp = xy_slurm_handler(config, sfiles)
} else {
  tmp = mclapply(sfiles, function(p){
    plot_xy(config, p)
  }, mc.cores = config$cores)
}


#summarize the spatial files
spatial_summary = mclapply(sfiles, function(p){
  generate_marker_summary(config, p)
}, mc.cores = config$cores) %>%
  do.call(bind_rows, .)
colnames(spatial_summary)[1] = config$variables$sample_id
fwrite(spatial_summary,
       file.path(config$paths$output, config$paths$sample))

#distribution plots
plot_distributions(config = config)

#calculate spatial metrics
cat("Calculating Metrics\n")
if(config$slurm){
  tmp = metric_slurm_handler(config, sfiles)
} else {
  tmp = mclapply(config$metrics, function(m){
    mclapply(sfiles, function(f){
      if(m == "kest") calculate_kest(config, f)
      if(m == "gest") calculate_gest(config, f)
      if(m == "dbscan") calculate_dbscan(config, f)
      if(m == "full_graph") full_interaction_graph(config, f)
    }, mc.allow.recursive = TRUE)
  }, mc.cores = config$cores, mc.allow.recursive = TRUE)
}

#plot functions for spatial metrics
cat("Plotting results\n")
cat("\tFinding derived results\n")
sm_files = lapply(config$metrics, function(m){
  list.files(file.path(config$paths$output, "metrics", m), full.names = TRUE)
})
names(sm_files) = config$metrics
cat("\tCreating plots\n")
if(config$slurm){
  tmp = mplot_slurm_handler(config, sm_files)
} else {
  tmp = mclapply(names(sm_files), function(m){
    tmp = mclapply(sm_files[[m]], function(f){
      if(m == "kest") plot_kest(config, f)
      if(m == "gest") plot_gest(config, f)
      if(m == "dbscan") plot_dbscan(config, f)#
      if(m == "full_graph") plot_full_graph(config, f)
    }, mc.allow.recursive = TRUE)
  }, mc.cores = config$cores, mc.allow.recursive = TRUE)
}



#doens't work on HPC - need jpeglib installed
#waiting on mehmet
# #paulos magic
# # combine all pdf files in /figures ---------------------------------------
# # Create a PDF file for the combined report
# pdf_output_path <- file.path(config$paths$output, "combined_report.pdf")
# output_dir = file.path(config$paths$output,
#                        "figures")
# 
# # Find all directories inside the output_dir
# dirs <- list.dirs(output_dir, full.names = TRUE, recursive = FALSE)
# 
# # Collect all PDF files in the directory, recursively
# pdf_files <- c()
# for (dir in dirs) {
#   # Find all PDF files in the current directory
#   pdf_files <- c(pdf_files, list.files(dir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE))
# }
# 
# # Check if there are any PDFs to combine
# if (length(pdf_files) > 0) {
#   cat("Found the following PDF files:\n")
#   print(pdf_files)
#   
#   # Combine the PDFs into one report using qpdf
#   qpdf::pdf_combine(pdf_files, pdf_output_path)
#   
#   cat("Saving combined PDF report to:", pdf_output_path, "\n")
# } else {
#   cat("No PDF files found to combine.\n")
# }
