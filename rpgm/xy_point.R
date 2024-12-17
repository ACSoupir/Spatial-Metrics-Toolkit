
plot_xy = function(config, path){
  #set correct path
  setwd(config$paths$root)
  #read in files
  df = fread(path)
  #convert data to long format
  id = unique(df[[config$variables$sample_id]])
  #get data
  df2 = df %>% 
    as.data.frame(check.names = FALSE) %>%
    select(any_of(c(config$variables$x_value, config$variables$y_value,
                    config$variables$markers,
                    config$variables$tissue_class_label))) %>%
    mutate(background = ifelse(rowSums(across(!!config$variables$markers)) == 0,
                               1, #not positive
                               0)) %>%
    pivot_longer(cols = c(!!config$variables$markers, background), 
                 values_to = "positive", names_to = "Marker") %>%
    filter(positive == 1)
  
  #plot
  if(!is.null(config$variables$tissue_class_label)){
    p = df2 %>%
      ggplot() + 
      geom_point(data = . %>%
                   filter(Marker == "background"),
                 aes(x = get(config$variables$x_value), y = get(config$variables$y_value), shape = get(config$variables$tissue_class_label)),
                 color = 'grey75') +
      geom_point(data = . %>%
                   filter(Marker != "background"),
                 aes(x = get(config$variables$x_value), y = get(config$variables$y_value), shape = get(config$variables$tissue_class_label), color = Marker)) +
      guides(shape=guide_legend(title=config$variables$tissue_class_label)) +
      coord_equal() +
      theme_classic()
  } else {
    p = df2 %>%
      ggplot() + 
      geom_point(data = . %>%
                   filter(Marker == "background"),
                 aes(x = get(config$variables$x_value), y = get(config$variables$y_value)),
                 color = 'grey75') +
      geom_point(data = . %>%
                   filter(Marker != "background"),
                 aes(x = get(config$variables$x_value), y = get(config$variables$y_value), color = Marker)) +
      coord_equal() +
      theme_classic()
  }
  
  pdf(file.path(config$paths$output, 'figures/point/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
      height = 10, width = 10)
  print(p +
          labs(x = config$variables$x_value,
               y = config$variables$y_value,
               title = id))
  dev.off()
  return("Success")
}