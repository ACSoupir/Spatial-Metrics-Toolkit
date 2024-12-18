#requires igraph
full_interaction_graph = function(config, path){
  #set correct path
  setwd(config$paths$root)
  #read in files
  df = fread(path)
  
  #markers
  markers = config$variables$markers %>% setNames(., .)
  radii = eval(parse(text = config$variables$radii_range))
  
  df_tmp = df %>%
    select(any_of(c(config$variables$x_value, 
                    config$variables$y_value,
                    markers, 
                    config$variables$tissue_class_label)))
  
  dist_mat = as.matrix(dist(df_tmp[,1:2]))
  
  # Identify interactions based on proximity threshold
  interactions <- which(dist_mat <= max(radii) & dist_mat > 0, arr.ind = TRUE)
  
  # Create a data frame of edges (interacting cells)
  edges <- data.table(
    source = interactions[, 1],
    target = interactions[, 2],
    distance = dist_mat[interactions]
  ) %>%
    mutate(x_source = df_tmp$x[source],
           y_source = df_tmp$y[source],
           x_target = df_tmp$x[target],
           y_target = df_tmp$y[target])
  #make metadata for the points
  res = list(
    Edges = edges,
    MetaData = df_tmp
  )
  saveRDS(res,
          file.path(config$paths$output, 'metrics/full_graph/', paste0(basename(gsub(".csv.*", "", path)), ".rds")))
  
  return("Success")
  
}

#dont think plotting this is really the best thing
plot_full_graph = function(config, path){
  #set correct path
  setwd(config$paths$root)
  #convert radii from character input to values
  radii = eval(parse(text = config$variables$radii_range))
  #identify a few to plot with an inverted probability (more low, less high)
  tmp_r_vec = log1p(seq(radii)) %>% round() #log2 the seq of radii vs range (ignores actual values)
  tmp_r_vec = radii[!duplicated(tmp_r_vec)]
  #read in data
  dat = readRDS(path)
  p_radii = tmp_r_vec[tmp_r_vec > min(dat$Edges$distance)]
  plots = lapply(p_radii, function(r){
    full_graph_ggplot(config, dat, radius = r)
  })
  pdf(file.path(config$paths$output, 'figures/metrics/full_graph/', paste0(basename(gsub(".rds.*", "", path)), ".pdf")),
      height = 7, width = 10)
  tmp = lapply(plots, print)
  dev.off()
  
  return("Success")
}

full_graph_ggplot = function(config, fig_dat, radius){ #full interaction graph
  #plot edges that fall under radius
  point_dat = fig_dat$MetaData %>%
    mutate(Other = ifelse(rowSums(across(!!config$variables$markers)) == 0,
                          1, #not positive
                          0)) %>%
    pivot_longer(cols = c(!!config$variables$markers, "Other"), 
                 values_to = "positive", names_to = "Marker") %>%
    filter(positive == 1)
  p1 = ggplot() + 
    geom_segment(data = fig_dat$Edges %>%
                   filter(distance < radius),
                 aes(x = x_source,
                     y = y_source,
                     xend = x_target,
                     yend = y_target), alpha = 0.3
                 ) 
  if(!is.null(config$variables$tissue_class_label)){
    p1 = p1 +
      geom_point(data = point_dat %>%
                          filter(Marker == "Other"),
                        aes(x = x, y = y,
                            shape = as.factor(get(config$variables$tissue_class_label))), 
                        color = 'grey65', alpha = 0.5) +
      geom_point(data = point_dat %>%
                   filter(Marker != "Other"),
                 aes(x = x, y = y, 
                     color = as.factor(Marker),
                     shape = as.factor(get(config$variables$tissue_class_label)))) +
      scale_shape_manual(values = c(16, 3), name = config$variables$tissue_class_label)
  } else {
    p1 = p1 +
      geom_point(data = point_dat %>%
                   filter(Marker == "Other"),
                 aes(x = x, y = y), 
                 color = 'grey65', alpha = 0.5) +
      geom_point(data = point_dat %>%
                   filter(Marker != "Other"),
                 aes(x = x, y = y, 
                     color = as.factor(Marker)))
  }
  p1 +
    coord_equal() +
    guides(color = guide_legend(title = "Marker")) +
    theme_classic() +
    labs(x = config$variables$x_value,
         y = config$variables$y_value,
         title = paste0("Radius = ", radius))
}
