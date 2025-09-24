
#NN G exact
#https://github.com/ACSoupir/PermFree-NearestNeighborG/
Rcpp::sourceCpp("src/compute_neighbor_counts_opt.cpp")
#faster matrix mulipliation
Rcpp::sourceCpp("src/mat_mult_eigen.cpp")

#combinatorial calculations in log space for large N samples
log_comb = function(n, k) {
  lgamma(n + 1) - (lgamma(k + 1) + lgamma(n - k + 1))
}

binom_prob = function(n, k, x) {
  exp(log_comb(n - x, k) - log_comb(n, k))
}


#calculate gest
calculate_gest_exactCSR = function(config, path){
  #set correct path
  setwd(config$paths$root)
  #read in files
  df = fread(path)
  #split based on the classifier label colun
  if(!is.null(config$variables$tissue_class_label)){
    df_list = split(df, df[[config$variables$tissue_class_label]])
  } else {
    df_list = list("Overall" = df)
  }
  #identify windows of compartments
  windows = lapply(df_list, function(x){
    convexhull.xy(x[[config$variables$x_value]], 
                  x[[config$variables$y_value]])
  })
  #point pattern objects
  pp_objs = lapply(seq(df_list), function(l){
    x = df_list[[l]]
    ppp(x[[config$variables$x_value]], 
        x[[config$variables$y_value]], window = windows[[l]],
        marks = x[,config$variables$markers, with = FALSE])
  })
  gest_res = lapply(seq(df_list), function(l){
    p = pp_objs[[l]]
    # x = df_list[[l]] %>%
    #   as.data.frame(check.names = FALSE)
    # id = unique(x[[config$variables$sample_id]])
    res = lapply(config$variables$markers, function(marker){
      #print(marker)
      p2 = subset(p, marks(p)[[marker]] == 1)
      if(npoints(p2) < 2){
        out = data.frame(r = eval(parse(text = config$variables$radii_range))) %>%
          mutate(Marker = marker, 
                 .before = 1)
      } else {
        out = Gest(p2, r=eval(parse(text = config$variables$radii_range)), correction = "none") %>% #don't know whether edge correction alters the results of if it cancels out
          as.data.frame(check.names = FALSE) %>%
          mutate(Marker = marker, 
                 .before = 1)
        pcells = npoints(p2)
        n_points = npoints(p)
        counts = compute_neighbor_counts_opt(p$x, p$y, eval(parse(text = config$variables$radii_range)))
        binom_vals = binom_prob(n_points-1, pcells-1, counts)
        summation2 = colSums(binom_vals)
        
        exact_csr_values = 1 - (1 / n_points) * summation2
        
        out$exactCSR = exact_csr_values
      }
      return(out)
    }) %>%
      do.call(bind_rows, .)
    
    #compartment if needed
    if(!is.null(config$variables$tissue_class_label)){
      res %>%
        mutate(!!config$variables$tissue_class_label := names(df_list)[l]) %>%
        return()
    } else {
      return(res)
    }
    
  }) %>%
    do.call(bind_rows, .) %>%
    dplyr::mutate(!!config$variables$sample_id := unique(df[[config$variables$sample_id]]),
                  .before = 1)
  fwrite(gest_res,
         file.path(config$paths$output, 'metrics/gest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".csv.gz")),
         compress = "gzip")
  
  return("Success")
}

#plot it
plot_gest_exactCSR = function(config, path){
  #set correct path
  setwd(config$paths$root)
  df = fread(path, data.table = FALSE) %>%
    mutate(r = as.numeric(r))
  if(!is.null(config$variables$tissue_class_label)){
    p = df %>%
      ggplot() +
      geom_line(aes(x = r, y = raw - exactCSR, color = Marker)) +
      facet_grid(get(config$variables$tissue_class_label)~.) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/gest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = 10)
    print(p)
    dev.off()
  } else {
    p = df %>%
      ggplot() +
      geom_line(aes(x = r, y = raw - exactCSR, color = Marker)) +
      theme_classic()
    pdf(file.path(config$paths$output, 'figures/metrics/gest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 6, width = 10)
    print(p)
    dev.off()
  }
  
  return("Success")
}
