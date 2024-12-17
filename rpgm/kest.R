
#calculate kest
calculate_kest = function(yaml, path){
  #read in files
  df = fread(path)
  #split based on the classifier label colun
  if(!is.null(yaml$variables$tissue_class_label)){
    df_list = split(df, df[[yaml$variables$tissue_class_label]])
  } else {
    df_list = list("Overall" = df)
  }
  #identify windows of compartments
  windows = lapply(df_list, function(x){
    convexhull.xy(x[[yaml$variables$x_value]], 
                  x[[yaml$variables$y_value]])
  })
  #point pattern objects
  pp_objs = lapply(seq(df_list), function(l){
    x = df_list[[l]]
    ppp(x[[yaml$variables$x_value]], 
        x[[yaml$variables$y_value]], window = windows[[l]],
        marks = x[,yaml$variables$markers, with = FALSE])
  })
  kest_res = lapply(seq(df_list), function(l){
    p = pp_objs[[l]]
    x = df_list[[l]] %>%
      as.data.frame(check.names = FALSE)
    id = unique(x[[yaml$variables$sample_id]])
    res = lapply(yaml$variables$markers, function(marker){
      #print(marker)
      p2 = subset(p, marks(p)[[marker]] == 1)
      if(npoints(p2) < 2){
        out = data.frame(r = eval(parse(text = yaml$variables$radii_range))) %>%
          mutate(!!yaml$variables$sample_id := id,
                 Marker = marker, 
                 .before = 1)
      } else {
        out = Kest(p2, r=eval(parse(text = yaml$variables$radii_range))) %>%
          as.data.frame(check.names = FALSE) %>%
          mutate(!!yaml$variables$sample_id := id,
                 Marker = marker, 
                 .before = 1)
      }
      
      if(!is.null(yaml$variables$tissue_class_label)){
        out %>%
          mutate(!!yaml$variables$tissue_class_label := names(df_list)[l]) %>%
          return()
      } else {
        return(out)
      }
     
    }) %>%
      do.call(bind_rows, .)
  }) %>%
    do.call(bind_rows, .)
  fwrite(kest_res,
         file.path(yaml$paths$output, 'metrics/kest/', paste0(basename(gsub(".csv.*", "", path)), ".csv.gz")),
         compress = "gzip")
  
  return("Success")
  
}

#plot it
plot_kest = function(config, path){
  df = fread(path, data.table = FALSE) %>%
    mutate(r = as.numeric(r))
  if(!is.null(config$variables$tissue_class_label)){
    p = df %>%
      ggplot() +
      geom_line(aes(x = r, y = border - theo, color = Marker)) +
      facet_grid(get(config$variables$tissue_class_label)~.) +
      theme_classic()
    
    pdf(file.path(config$paths$output, 'figures/metrics/kest/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = 10)
    print(p)
    dev.off()
  } else {
    p = df %>%
      ggplot() +
      geom_line(aes(x = r, y = border - theo, color = Marker)) +
      theme_classic()
    pdf(file.path(config$paths$output, 'figures/metrics/kest/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 6, width = 10)
    print(p)
    dev.off()
  }
  
  return("Success")
}
