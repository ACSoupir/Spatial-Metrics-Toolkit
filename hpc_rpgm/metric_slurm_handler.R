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
                       "calculate_kest_exactCSR",
                       "calculate_gest",
                       "calculate_gest_exactCSR",
                       "calculate_dbscan",
                       "full_interaction_graph",
                       "calculate_kcross",
                       "calculate_kcross_exact",
                       "calculate_gcross",
                       "calculate_gcross_exact")
  )
  
  res = get_slurm_out(metric_job, outtype = "raw")
  res_files = list.files(".", pattern = "results") %>%
    grep("RDS", ., value = TRUE)
  file.remove(res_files)
  return("Done")
}

metric_slurm_shuffler = function(metric, config, path){
  if(metric == "kest") out = calculate_kest(config, path)
  if(metric == "kest_exactCSR") out = calculate_kest_exactCSR(config, path)
  if(metric == "gest") out = calculate_gest(config, path)
  if(metric == "gest_exactCSR") out = calculate_gest_exactCSR(config, path)
  if(metric == "dbscan") out = calculate_dbscan(config, path)
  if(metric == "full_graph") out = full_interaction_graph(config, path)
  
  #bivariate
  if(metric == "kcross") out = calculate_kcross(config, path)
  if(metric == "kcross_exactCSR") out = calculate_kcross_exact(config, path)
  if(metric == "gcross") out = calculate_gcross(config, path)
  if(metric == "gcross_exactCSR") out = calculate_gcross_exact(config, path)
  return(out)
}