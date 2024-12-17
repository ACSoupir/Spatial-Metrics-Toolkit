#calculate the metric
calculate_dbscan = function(config, path){
  #set correct path
  setwd(config$paths$root)
  #read in files
  df = fread(path, data.table = FALSE)
  #split based on the classifier label colun
  if(!is.null(config$variables$tissue_class_label)){
    df_list = split(df, df[[config$variables$tissue_class_label]])
  } else {
    df_list = list("Overall" = df)
  }
  
  #broken down by
  #tissue_class_label > marker > radius
  radii = eval(parse(text = config$variables$radii_range)) %>%
    setNames(.,.)
  out = lapply(df_list, function(tcl){
    res = lapply(config$variables$markers, function(marker){
      #print(marker)
      #subset to positive for marker
      wkn_df = tcl %>% filter(get(marker) == 1)
      #handle 0
      if(nrow(wkn_df) < 2) return() #do something to skip/handle
      x = as.matrix(wkn_df[,c(config$variables$x_value,
                              config$variables$y_value)])
      db_out = lapply(radii, function(r){
        set.seed(333)
        db = dbscan(x, eps = r, minPts = 5)
        cluster_sizes = table(db$cluster)
        cluster_props = cluster_sizes / sum(cluster_sizes)
        return(list(clusters = db$cluster,
                    num_clusters = length(unique(db$cluster)) - 1,
                    cluster_sizes = cluster_sizes,
                    cluster_props = cluster_props))
      })
    })
    names(res) = config$variables$markers
    return(res)
  })
  saveRDS(out,
          file.path(config$paths$output, 'metrics/dbscan/', paste0(basename(gsub(".csv.*", "", path)), ".rds")))
  
  return("Success")
}

#plot it
plot_dbscan = function(config, path){
  #set correct path
  setwd(config$paths$root)
  db_res = readRDS(path)
  #convert results to table
  out = lapply(names(db_res), function(t){
    tis = db_res[[t]]
    tmp = lapply(names(tis), function(m){
      marker = tis[[m]]
      lapply(marker, function(r){
        data.table(clusters = r$num_clusters)
      }) %>% do.call(bind_rows, .) %>% 
        rownames_to_column("radius") %>%
        mutate(Marker = m)
    }) %>% do.call(bind_rows, .)
    if(!is.null(config$variables$tissue_class_label)){
      tmp %>%
        mutate(!!config$variables$tissue_class_label := t) %>%
        return()
    } else {
      return(tmp)
    }
  }) %>%
    do.call(bind_rows, .)
  #create plot
  p = out %>%
    mutate(radius = as.numeric(radius)) %>%
    ggplot() + 
    geom_line(aes(x = radius, y = clusters, color = Marker)) +
    theme_classic()
  
  if(!is.null(config$variables$tissue_class_label)){
    p = p +
      facet_grid(get(config$variables$tissue_class_label)~.)
    
    pdf(file.path(config$paths$output, 'figures/metrics/dbscan/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = 10)
    print(p)
    dev.off()
  } else {
    pdf(file.path(config$paths$output, 'figures/metrics/dbscan/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
        height = 7, width = 10)
    print(p)
    dev.off()
  }
  return("Success")
}