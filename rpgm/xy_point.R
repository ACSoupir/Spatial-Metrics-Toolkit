
plot_xy = function(yaml, path){
  #read in files
  df = data.table::fread(path)
  #convert data to long format
  id = unique(df[[yaml$variables$sample_id]])
  #get data
  df2 = df %>% 
    as.data.frame(check.names = FALSE) %>%
    select(any_of(c("x", "y",
                    yaml$variables$markers,
                    yaml$variables$tissue_class_label))) %>%
    mutate(background = ifelse(rowSums(across(!!yaml$variables$markers)) == 0,
                               1, #not positive
                               0)) %>%
    pivot_longer(cols = c(!!yaml$variables$markers, background), 
                 values_to = "positive", names_to = "Marker") %>%
    filter(positive == 1)
  
  #plot
  p = df2 %>%
    ggplot() + 
    geom_point(data = . %>%
                 filter(Marker == "background"),
               aes(x = x, y = y, shape = get(yaml$variables$tissue_class_label)),
               color = 'grey75') +
    geom_point(data = . %>%
                 filter(Marker != "background"),
               aes(x = x, y = y, shape = get(yaml$variables$tissue_class_label), color = Marker)) +
    guides(shape=guide_legend(title=yaml$variables$tissue_class_label)) +
    coord_equal() +
    theme_classic()
  
  suppressMessages(ggsave(file.path(yaml$paths$output, 'figures/point', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
         height = 10, width = 10))
  # pdf(file.path(yaml$paths$output, 'figures/point/', paste0(basename(gsub(".csv.*", "", path)), ".pdf")),
  #     height = 10, width = 10)
  # print(p)
  # dev.off()
}