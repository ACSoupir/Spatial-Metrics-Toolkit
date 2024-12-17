#This function will generated the plot of cell type distribution across samples
# based on the mif summary file. The interested cell types were setup at yaml file  

#libraries needed
# require(ggplot2)
# require(Polychrome) # generated colors
# require(ggpubr)
# cell_distribution_plot <- function(yaml_config=NULL){
#   # functions for barplot
#   box_plots <- function(data =NULL, 
#                         markers = NULL, 
#                         label =NULL){
#     plot_data <- matrix(0,nrow=nrow(data) * length(markers),ncol=3)
#     m <- 0
#     for(i in markers){
#       for (j in 1:nrow(data)){
#         m <- m + 1
#         plot_data[m,1] <- data[j,1] # sample name
#         plot_data[m,2] <- data[j,i] # Percentage of positive cells
#         plot_data[m,3] <- i # cell types
#       }
#     }
#     plot_data <- data.frame(plot_data)
#     plot_data[,1] <- gsub(".tif", "", plot_data[,1])
#     plot_data[,3] <- gsub("^X..", "",plot_data[,3])
#     plot_data[,3] <- gsub(".Cells", "",plot_data[,3])
#     colors <- createPalette(length(unique(plot_data[,1])), c("#2a6ebb", "#de3831", "#007367"), range = c(30, 80))
#     names(colors) <- unique(plot_data[,1])
#     theme_set(theme_bw())
#     ggplot(plot_data, aes(y=as.numeric(plot_data[,2]), x= plot_data[,3])) +
#       geom_boxplot(show.legend = F) +
#       geom_point(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.5),
#                  aes(color = factor(plot_data[,1])), show.legend =T) +
#       labs(x = '', 
#            y = "Percentage of positive cells (%)",
#            color = "Sample ID"
#       ) +
#       scale_x_discrete(guide = guide_axis(angle = 45)) +
#       scale_color_manual(values=colors) +
#       theme(
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(face='bold', size=10),
#         axis.text.y = element_text(face='bold', size=10),
#         axis.title.y = element_text(face='bold', size=10),
#         axis.title.x = element_text(face='bold', size=10),
#         legend.title = element_text(face='bold'),
#         legend.text = element_text(face='bold'))
#     
#     if (is.null(label)){
#     out1 <- paste0(outputFilePath, "/", "cell_distribution_acorss_sample_boxplot.pdf", sep="")
#     }else{
#       out1 <- paste0(outputFilePath, "/", label, "_cell_distribution_acorss_sample_boxplot.pdf", sep="")
#     }
#     ggsave(out1)
#     }
#   # summary file data read and generated the figures 
#   if(is.null(yaml_config$variables$sample_id) & 
#      is.null(yaml_config$variables$markers) & 
#      is.null(yaml_config$paths$output) &
#      is.null(yaml_config$paths$sample) |
#      is.null(yaml_config$variables$tissue_class_label)){
#     stop("sample_id, marker set, summary file path, and output path 
#     must be provided or tissue_class_label in optional in yaml config file.", call.=FALSE)}
#   summary_path <- file.path(yaml_config$paths$output, yaml_config$paths$sample)
#   marker_set <- yaml_config$variables$markers
#   sample_col <- yaml_config$variables$sample_id
#   outputFilePath = yaml_config$paths$output
#   data <- read.csv(file = summary_path, header = TRUE)
#   # marker_set1 <- gsub(".Positive", "", marker_set)
# #  marker_set1 <- paste0("X..", marker_set1, ".Positive.Cells")
#   marker_set1 <- paste0("X..", marker_set, ".Cells")
#   if (is.null(yaml_config$variables$tissue_class_label)){
#     data1 <- data[, c(sample_col, marker_set1)]
#     box_plots(data1, marker_set1, label = NULL)
#   } else {
#     class_label <- yaml_config$variables$tissue_class_label
#    for (i in unique(data[,class_label])){
#      data1 <- data[data[,class_label]== i, c(sample_col, marker_set1)]
#      box_plots(data1, marker_set1, label = i)
#    } 
#   }
# }

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
