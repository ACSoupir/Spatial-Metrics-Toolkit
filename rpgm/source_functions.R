
# sourcing functions ------------------------------------------------------

#processing functions
source("rpgm/create_folders.R")
source("rpgm/kest.R")
source("rpgm/gest.R")
source("rpgm/dbscan.R")
source("rpgm/xy_point.R")
source("rpgm/Summarise_function.R")
source("rpgm/validations.R")
source("rpgm/0.4.cell_distribution_across_samples_boxplot.R")

#slurm helper functions
source("hpc_rpgm/get_slurm_out.R")
source("hpc_rpgm/xy_slurm_handler.R")
source("hpc_rpgm/metric_slurm_handler.R")
source("hpc_rpgm/m-plot_slurm_handler.R")
