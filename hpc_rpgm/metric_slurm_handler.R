metric_slurm_handler = function(config, sfiles){
  params = expand_grid(metric = config$metrics,
                       config = I(list(config)),
                       path = sfiles)
  
  metric_job = slurm_apply(
    metric_slurm_shuffler,
    params,
    jobname = "rslurm_metric",
    nodes = nrow(params),
    cpus_per_node = 1,
    global_objects = c("calculate_kest",
                       "calculate_gest",
                       "calculate_dbscan",
                       "full_interaction_graph")
  )
  
  res = get_slurm_out(metric_job, outtype = "raw")
  res_files = list.files(".", pattern = "results") %>%
    grep("RDS", ., value = TRUE)
  file.remove(res_files)
  return("Done")
}

metric_slurm_shuffler = function(metric, config, path){
  if(metric == "kest") out = calculate_kest(config, path)
  if(metric == "gest") out = calculate_gest(config, path)
  if(metric == "dbscan") out = calculate_dbscan(config, path)
  if(metric == "full_graph") out = full_interaction_graph(config, path)
  return(out)
}