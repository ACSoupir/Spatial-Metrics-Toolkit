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
                       "plot_kest_exactCSR",
                       "plot_gest",
                       "plot_gest_exactCSR",
                       "plot_dbscan",
                       "plot_full_graph", "full_graph_ggplot",
                       "plot_kcross", "plot_kcross_exact",
                       "plot_gcross", "plot_gcross_exact")
  )
  
  res = get_slurm_out(mplot_job, outtype = "raw")
  res_files = list.files(".", pattern = "results") %>%
    grep("RDS", ., value = TRUE)
  file.remove(res_files)
  return("Done")
}

mplot_slurm_shuffler = function(metric, config, path){
  if(metric == "kest") out = plot_kest(config, path)
  if(metric == "kest_exactCSR") out = plot_kest_exactCSR(config, path)
  if(metric == "gest") out = plot_gest(config, path)
  if(metric == "gest_exactCSR") out = plot_gest_exactCSR(config, path)
  if(metric == "dbscan") out = plot_dbscan(config, path)
  if(metric == "full_graph") out = plot_full_graph(config, path)
  #bivariate
  if(metric == "kcross") out = plot_kcross(config, path)
  if(metric == "kcross_exactCSR") out = plot_kcross_exact(config, path)
  if(metric == 'gcross') out = plot_gcross(config, path)
  if(metric == "gcross_exactCSR") out = plot_gcross_exact(config, path)
  return(out)
}