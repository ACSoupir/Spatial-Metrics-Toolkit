mplot_slurm_handler = function(config, sm_files){
  params = as.data.frame(sm_files) %>% 
    pivot_longer(!!names(.), names_to = "metric", values_to = "path") %>%
    mutate(config = I(list(config)))
  
  mplot_job = slurm_apply(
    mplot_slurm_shuffler,
    params,
    jobname = "rslurm_mplot",
    nodes = nrow(params),
    cpus_per_node = 1,
    global_objects = c("plot_kest",
                       "plot_gest",
                       "plot_dbscan")
  )
  
  res = get_slurm_out(mplot_job, outtype = "raw")
  res_files = list.files(".", pattern = "results") %>%
    grep("RDS", ., value = TRUE)
  file.remove(res_files)
  return("Done")
}

mplot_slurm_shuffler = function(metric, config, path){
  if(metric == "kest") out = plot_kest(config, path)
  if(metric == "gest") out = plot_gest(config, path)
  if(metric == "dbscan") out = plot_dbscan(config, path)
  return(out)
}