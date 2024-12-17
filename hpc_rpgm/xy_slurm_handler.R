xy_slurm_handler = function(config, sfiles){
  params = expand_grid(config = I(list(config)),
                       path = sfiles)
  
  xy_job = slurm_apply(
    plot_xy,
    params,
    jobname = "rslurm_xy",
    nodes = nrow(params),
    cpus_per_node = 1
  )
  
  res = get_slurm_out(xy_job, outtype = "raw")
  res_files = list.files(".", pattern = "results") %>%
    grep("RDS", ., value = TRUE)
  file.remove(res_files)
  return("Done")
}