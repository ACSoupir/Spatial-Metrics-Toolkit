
plot_distributions = function(config){
  df = fread(file.path(config$paths$output, config$paths$sample), data.table = FALSE) %>%
    select(any_of(config$variables$sample_id),
           any_of(c(contains('%'))),
           any_of(config$variables$tissue_class_label)) %>%
    pivot_longer(cols = contains(config$variables$markers),
                 names_to = "Marker Abundance", values_to = "Percent")
  
  p = df %>%
    ggplot() + 
    geom_boxplot(aes(x = `Marker Abundance`, y = Percent)) +
    geom_jitter(aes(x = `Marker Abundance`, y = Percent, color = get(config$variables$sample_id))) +
    guides(color=guide_legend(title=config$variables$sample_id)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if(!is.null(config$variables$tissue_class_label)){
    p = p +
      facet_grid(get(config$variables$tissue_class_label)~.)
    pdf(file.path(config$paths$output, 'figures/barplot/marker_distribution.pdf'),
        height = 7, width = 10)
    print(p)
    dev.off()
  } else {
    pdf(file.path(config$paths$output, 'figures/barplot/marker_distribution.pdf'),
        height = 5, width = 10)
    print(p)
    dev.off()
  }
}
