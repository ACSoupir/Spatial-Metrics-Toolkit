

#' Title: Generate a graph representation of the spatial data 
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#' @param radius the defined radius for connecting the cells
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generateFullGraphPlotAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generateFullGraphPlotAsPDF <- function(yamlfile, csvfilePath, radius = 20) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  shape = categories

  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/InteractiveNetork_", clean_name, "_", marker, ".pdf")
    outputGraphname = paste0(outputFilePath, "/InteractiveNetork_", clean_name, "_", marker, ".graphml")
    
    data_base[which(data_base[,marker] == 1), marker] = "pos"
    data_base[which(data_base[,marker] == 0), marker] = "neg"
    
    plot_data = cbind(data_base$x, data_base$y, paste0(data_base[,marker], "_", data_base[,shape]))
    maker_shape =  paste0(marker, "_", shape)
    colnames(plot_data) = c("X", "Y", "cell_type")
    plot_data = as.data.frame(plot_data)
    plot_data[1:2] <- lapply(plot_data[1:2], as.numeric) # forcing the x and y as numeric
    
    mif_data <- plot_data
    mif_data$Cell_ID = rownames(mif_data)
    
    # Set distance threshold for proximity (e.g., 20 microns)
    proximity_threshold <- radius
    
    # Calculate pairwise distances
    
    distance_matrix <- as.matrix(dist(mif_data[, c("X", "Y")]))
    
    # Identify interactions based on proximity threshold
    interactions <- which(distance_matrix <= proximity_threshold & distance_matrix > 0, arr.ind = TRUE)
    
    # Create a data frame of edges (interacting cells)
    edges <- data.table(
      source = mif_data$Cell_ID[interactions[, 1]],
      target = mif_data$Cell_ID[interactions[, 2]],
      distance = distance_matrix[interactions]
    )
    
    # Generate an interaction network using igraph
    interaction_network <- graph_from_data_frame(edges, directed = FALSE)
    
    # Add cell type attributes to the network
    V(interaction_network)$cell_type <- mif_data$cell_type[match(V(interaction_network)$name, mif_data$Cell_ID)]
    
    # Plot the interaction network
    #plot(interaction_network,
    #     vertex.color = as.factor(V(interaction_network)$CD3..FOXP3._Classifier.Label),
    #     vertex.size = 2,
    #     vertex.label = NA,
    #     edge.width = 0.5,
    #     main = "Cell-Cell Proximity Interaction Network")
    
    # Optional: Export the network as a file
    write_graph(interaction_network, file = outputGraphname, format = "graphml")
    
    # Visualization: ggplot scatter plot with interactions
    p = ggplot(mif_data, aes(X, Y)) +
      geom_point(aes(color = cell_type), size = 2) +
      geom_segment(data = edges,
                   aes(x = mif_data$X[match(source, mif_data$Cell_ID)],
                       y = mif_data$Y[match(source, mif_data$Cell_ID)],
                       xend = mif_data$X[match(target, mif_data$Cell_ID)],
                       yend = mif_data$Y[match(target, mif_data$Cell_ID)]),
                   color = "gray", alpha = 0.5) +
      labs(title = "Cell-Cell Proximity Interaction Map",
           x = "X Coordinate", y = "Y Coordinate") +
      theme_minimal()
    
    ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}
