
#calculate kest
calculate_kest_exactCSR = function(config, path){
  #set correct path
  if(!(getwd() == config$paths$root)){
    setwd(config$paths$root)
  }
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
    marks = x[,config$variables$markers, with = FALSE]
    marks[] = lapply(marks, as.factor)
    ppp(x[[config$variables$x_value]], 
        x[[config$variables$y_value]], window = windows[[l]],
        marks = marks)
  })
  #calucalte k with exact CSR
  kest_res = lapply(seq(df_list), function(l){
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
        obs = Kest(p2, r=eval(parse(text = config$variables$radii_range))) %>%
          as.data.frame(check.names = FALSE) %>%
          select(-theo) %>%
          mutate(Marker = marker, 
                 .before = 1)
        exact = Kest(p, r=eval(parse(text = config$variables$radii_range))) %>%
          as.data.frame(check.names = FALSE) %>%
          select(-theo) %>%
          rename_at(-1, ~ paste0(.x, "_exact")) %>%
          mutate(Marker = marker, 
                 .before = 1)
        out = full_join(obs, exact, by = join_by(Marker, r))
      }
     
    }) %>%
      do.call(bind_rows, .)
    
    
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
  fwrite(kest_res,
         file.path(config$paths$output, 'metrics/kest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".csv.gz")),
         compress = "gzip")
  
  return("Success")
  
}

#plot it
plot_kest_exactCSR = function(config, path){
  #set correct path
  setwd(config$paths$root)
  df = fread(path, data.table = FALSE) %>%
    mutate(r = as.numeric(r))
  observed_cols = c('border', 'trans', 'iso')
  observed_cols = observed_cols[observed_cols %in% colnames(df)]
  for(cname in observed_cols){
    df[[paste0(cname, "_adjusted")]] = df[[cname]] - df[[paste0(cname, "_exact")]]
  }
  if(!is.null(config$variables$tissue_class_label)){
    p = df %>%
      pivot_longer(cols = contains("adjusted"), names_to = "edge_correction", values_to = "clustering") %>%
      ggplot() + 
      geom_line(aes(x = r, y = clustering, color = Marker, linetype = edge_correction)) +
      facet_grid(get(config$variables$tissue_class_label)~.) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/kest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = 10)
    print(p)
    dev.off()
  } else {
    p = df %>%
      pivot_longer(cols = contains("adjusted"), names_to = "edge_correction", values_to = "clustering") %>%
      ggplot() + 
      geom_line(aes(x = r, y = clustering, color = Marker, linetype = edge_correction)) +
      theme_classic()
    pdf(file.path(config$paths$output, 'figures/metrics/kest_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 6, width = 10)
    print(p)
    dev.off()
  }
  
  return("Success")
}
