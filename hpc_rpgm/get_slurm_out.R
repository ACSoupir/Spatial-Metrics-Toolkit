get_wd = function(wd = "NULL"){
  setwd(wd)
  
  v = system("pwd", intern = TRUE)
  return(v)
}

get_slurm_out = function (slr_job, outtype = "raw", wait = TRUE, ncores = NULL) 
{
  if (!(inherits(slr_job, "slurm_job"))) {
    stop("slr_job must be a slurm_job")
  }
  outtypes <- c("table", "raw")
  if (!(outtype %in% outtypes)) {
    stop(paste("outtype should be one of:", paste(outtypes, 
                                                  collapse = ", ")))
  }
  if (!(is.null(ncores) || (is.numeric(ncores) && length(ncores) == 
                            1))) {
    stop("ncores must be an integer number of cores")
  }
  if (wait) {
    rslurm:::wait_for_job(slr_job)
  }
  res_files <- paste0("results_", 0:(slr_job$nodes - 1), ".RDS")
  tmpdir <- "."
  missing_files <- setdiff(res_files, dir(path = tmpdir))
  if (length(missing_files) > 0) {
    missing_list <- paste(missing_files, collapse = ", ")
    warning(paste("The following files are missing:", missing_list))
  }
  res_files <- file.path(tmpdir, setdiff(res_files, missing_files))
  if (length(res_files) == 0) 
    return(NA)
  if (is.null(ncores)) {
    slurm_out <- lapply(res_files, readRDS)
  }
  else {
    slurm_out <- mclapply(res_files, readRDS, mc.cores = ncores)
  }
  slurm_out <- do.call(c, slurm_out)
  if (outtype == "table") {
    slurm_out <- as.data.frame(do.call(rbind, slurm_out))
  }
  slurm_out
}


# sjob1 <- slurm_apply(
#   get_wd,
#   params,
#   jobname="rslurm-test",
#   nodes=1,
#   cpus_per_node=1
# )
# 
# res <- get_slurm_out(sjob1, outtype="raw", wait = TRUE)
# 
# sjob1 = list()
# sjob1$jobname = "rslurmtest"
# sjob1$jobid = "11979848"
# sjob1$nodes = 1
# sjob1 = as(sjob1, "slurm_job")
