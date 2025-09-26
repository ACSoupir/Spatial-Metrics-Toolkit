
#calculate kest
calculate_kcross_exact = function(config, path){
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
    marks = x[,unique(unlist(config$variables$bivar_pairs)), with = FALSE]
    marks[] = lapply(marks, as.factor) #really is unneeded but YOLO
    ppp(x[[config$variables$x_value]], 
        x[[config$variables$y_value]], window = windows[[l]],
        marks = marks)
  })
  kest_res = lapply(seq(df_list), function(l){
    p = pp_objs[[l]]
    exact = Kest(p, r=eval(parse(text = config$variables$radii_range)), nlarge = 10000) %>% #forced 10k top limit
      as.data.frame(check.names = FALSE) %>%
      select(-theo) %>%
      rename_at(-1, ~ paste0(.x, "_exact"))
    #pairs
    res = lapply(config$variables$bivar_pairs, function(markers){
      #find cells to subset to
      mark_cells = marks(p) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), ~as.numeric(as.character(.x))), #forcing 1/0 for positive negative
                      id = 1:dplyr::n()) %>%
        dplyr::filter(get(markers[1]) != get(markers[2])) %>%
        dplyr::select(!!markers, id) #selects in correct order matching markers - anchor, counted
      counts = colSums(mark_cells[,markers, with = FALSE])
      if(any(counts < 2)){
        out = data.frame(r = eval(parse(text = config$variables$radii_range)))  %>%
          mutate(Anchor = markers[1],
                 Counted = markers[2], 
                 .before = 1)
      } else {
        cells = mark_cells %>%
          tidyr::pivot_longer(dplyr::all_of(markers), names_to = "marker", values_to = "positive") %>%
          dplyr::filter(positive == 1)
        p2 = subset(p, cells$id)
        marks(p2) = factor(cells$marker, levels = markers)
        obs = Kcross(p2, r=eval(parse(text = config$variables$radii_range))) %>%
          as.data.frame(check.names = FALSE) %>%
          mutate(Anchor = markers[1],
                 Counted = markers[2], 
                 .before = 1)
        out = full_join(obs, exact, join_by(r))
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
         file.path(config$paths$output, 'metrics/kcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".csv.gz")),
         compress = "gzip")
  
  return("Success")
  
}

#plot it
plot_kcross_exact = function(config, path){
  #set correct path
  setwd(config$paths$root)
  df = fread(path, data.table = FALSE) %>%
    mutate(r = as.numeric(r))
  observed_cols = c('border', 'trans', 'iso')
  observed_cols = observed_cols[observed_cols %in% colnames(df)]
  for(cname in observed_cols){
    df[[paste0(cname, "_adjusted")]] = df[[cname]] - df[[paste0(cname, "_exact")]]
  }
  c_num = length(unique(df$Anchor))
  if(!is.null(config$variables$tissue_class_label)){
    p = df %>%
      pivot_longer(cols = contains("adjusted"), names_to = "edge_correction", values_to = "clustering") %>%
      ggplot() + 
      geom_line(aes(x = r, y = clustering, color = Counted, linetype = edge_correction)) + 
      facet_grid(get(config$variables$tissue_class_label) ~ Anchor) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/kcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = (5*c_num)+0.5)
    print(p)
    dev.off()
  } else {
    c_num = length(unique(df$Anchor))
    p = df %>%
      pivot_longer(cols = contains("adjusted"), names_to = "edge_correction", values_to = "clustering") %>%
      ggplot() + 
      geom_line(aes(x = r, y = clustering, color = Counted, linetype = edge_correction)) + 
      facet_grid(. ~ Anchor) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/kcross_exactCSR/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = (5*c_num)+0.5)
    print(p)
    dev.off()
  }
  
  return("Success")
}
