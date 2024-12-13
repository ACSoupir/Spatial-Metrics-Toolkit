#run with
#"Rscript main.R -y config.yml"

#global libraries
cat("Prepping environment\n")
# List of packages
packages <- c("yaml", "optparse", "parallel", "data.table",
              "ggplot2", "tibble", "dplyr", "tidyr", "Polychrome", "ggpubr",
              "spatstat.geom", "spatstat.explore", "dbscan")

# Install missing packages 
missing_packages <- setdiff(packages,names(installed.packages()[,1]))
install.packages(missing_packages, quiet = TRUE)

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
    help = "number of cores for parallel processing",
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

#source validation functions
source("rpgm/validations.R")
cat("Running validations\n")
source("rpgm/validations.R")
validate_config(config)
# example call: function to validate df columns
# validate_df_columns(df, c("col1", "col2"))

#spatial files
cat("Identifying spatial files\n")
sfiles = list.files(config$paths$spatial, 
                    pattern = 'csv',
                    full.names = TRUE)

#source functions
cat("Sourcing functions\n")
source("rpgm/create_folders.R")
source("rpgm/kest.R")
source("rpgm/gest.R")
source("rpgm/dbscan.R")
source("rpgm/xy_point.R")
source("rpgm/Summarise_function.R")


#prep folder structure
cat("Creating folder structure\n")
createfolders(config)

#plot functions for raw csv inputs
#use data.table to import gz
cat("Plotting xy points\n")
tmp = mclapply(sfiles, function(p){
  plot_xy(config, p)
}, mc.cores = opt$cores)

#summarize the spatial files
spatial_summary = mclapply(sfiles, function(p){
  generate_marker_summary(config, p)
}, mc.cores = opt$cores) %>%
  do.call(bind_rows, .)
colnames(spatial_summary)[1] = config$variables$sample_id
fwrite(spatial_summary,
       file.path(config$paths$output, config$paths$sample))

#distribution plots
plot_distributions(config = config)

#calculate spatial metrics
cat("Calculating Metrics\n")
tmp = mclapply(config$metrics, function(m){
  mclapply(sfiles, function(f){
    if(m == "kest") calculate_kest(config, f)
    if(m == "gest") calculate_gest(config, f)
    if(m == "dbscan") calculate_dbscan(config, f)
  }, mc.allow.recursive = TRUE)
}, mc.cores = opt$cores, mc.allow.recursive = TRUE)

#plot functions for spatial metrics
cat("Plotting results\n")
cat("\tFinding derived results\n")
sm_files = lapply(config$metrics, function(m){
  list.files(file.path(config$paths$output, "metrics", m), full.names = TRUE)
})
names(sm_files) = config$metrics
cat("\tCreating plots\n")
tmp = mclapply(names(sm_files), function(m){
  tmp = mclapply(sm_files[[m]], function(f){
    if(m == "kest") plot_kest(config, f)
    if(m == "gest") plot_gest(config, f)
    if(m == "dbscan") plot_dbscan(config, f)
  }, mc.allow.recursive = TRUE)
}, mc.cores = opt$cores, mc.allow.recursive = TRUE)

