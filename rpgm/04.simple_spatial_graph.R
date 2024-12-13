

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
generateFullGraphPlotAsPDF <- function(yamlfile, csvfilePath, radius = 30) {
  
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

#' Title: Generate a sub clique network graph representation of the spatial data 
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#' @param radius the defined radius for connecting the cells
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generateSubCliqueGraphPlotAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generateSubCliqueGraphPlotAsPDF <- function(yamlfile, csvfilePath, radius = 30) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  shape = categories
  
  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/SubCliqueNetwork_", clean_name, "_", marker, ".pdf")
    outputGraphname = paste0(outputFilePath, "/SubCliqueNetwork_", clean_name, "_", marker, ".graphml")
    
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
    
    # Find cliques in the interaction network
    cliques <- cliques(interaction_network)
    
    # Filter cliques by size (e.g., only cliques of size >= 3)
    large_cliques <- Filter(function(c) length(c) >= 3, cliques)
    
    # Extract nodes involved in large cliques
    clique_nodes <- unique(unlist(large_cliques))
    
    # Create a subgraph of only the nodes involved in cliques
    clique_subgraph <- induced_subgraph(interaction_network, vids = clique_nodes)
    
    write_graph(clique_subgraph, file = outputGraphname, format = "graphml")
    
    # Extract node and edge data from the clique_subgraph
    clique_nodes_df <- data.frame(
      name = V(clique_subgraph)$name,
      cell_type = V(clique_subgraph)$cell_type,
      X = mif_data$X[match(V(clique_subgraph)$name, mif_data$Cell_ID)],
      Y = mif_data$Y[match(V(clique_subgraph)$name, mif_data$Cell_ID)]
    )
    
    
    
    clique_edges_df <- as.data.frame(get.data.frame(clique_subgraph, what = "edges"))
    
    # Add coordinates for edges
    clique_edges_df <- clique_edges_df %>%
      mutate(
        x_start = clique_nodes_df$X[match(from, clique_nodes_df$name)],
        y_start = clique_nodes_df$Y[match(from, clique_nodes_df$name)],
        x_end = clique_nodes_df$X[match(to, clique_nodes_df$name)],
        y_end = clique_nodes_df$Y[match(to, clique_nodes_df$name)]
      )
    
    p = ggplot() +
      # Add edges
      geom_segment(data = clique_edges_df,
                   aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                   color = "gray", alpha = 0.7) +
      # Add nodes
      geom_point(data = clique_nodes_df,
                 aes(x = X, y = Y, color = cell_type),
                 size = 2) +
      labs(title = "Cliques in the Interaction Network",
           x = "X Coordinate", y = "Y Coordinate", color = "Cell Type") +
      theme_minimal()
    
     ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}



#' Title: Identify the central hubs based on the interaction network of the spatial data
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#' @param radius the defined radius for connecting the cells
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generateCentralityPlotAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generateCentralityPlotAsPDF <- function(yamlfile, csvfilePath, radius = 30) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  shape = categories
  
  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/CentralityNetwork_", clean_name, "_", marker, ".pdf")
    outputGraphname = paste0(outputFilePath, "/CentralityNetwork_", clean_name, "_", marker, ".graphml")
    
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
    
    # Calculate degree centrality
    degree_centrality <- degree(interaction_network)
    
    # Add centrality values to the node data
    node_data <- data.frame(
      name = V(interaction_network)$name,
      X = mif_data$X[match(V(interaction_network)$name, mif_data$Cell_ID)],
      Y = mif_data$Y[match(V(interaction_network)$name, mif_data$Cell_ID)],
      cell_type = mif_data$cell_type[match(V(interaction_network)$name, mif_data$Cell_ID)],
      degree_centrality = degree_centrality
    )
    
    
    # Plot centrality using ggplot2
    library(ggplot2)
    
    p = ggplot(node_data, aes(X, Y)) +
      geom_point(aes(size = degree_centrality, alpha = degree_centrality, color = cell_type), alpha = 0.7) +
      scale_size_continuous(name = "Degree Centrality") +
      scale_alpha_continuous(name = "Centrality Opacity", range = c(0.0, 0.3)) +
      labs(title = "Centrality Plot", x = "X Coordinate", y = "Y Coordinate") +
      theme_minimal()
    
    
    ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}


#' Title: Generate a betweenness measure of the proximity network of the spatial data
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#' @param radius the defined radius for connecting the cells
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generateBetweenessMeasurePlotAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generateBetweenessMeasurePlotAsPDF <- function(yamlfile, csvfilePath, radius = 30) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  shape = categories
  
  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/BetweennessNetwork_", clean_name, "_", marker, ".pdf")
    outputGraphname = paste0(outputFilePath, "/BetweenessNetwork_", clean_name, "_", marker, ".graphml")
    
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
    
    
    # Calculate betweenness centrality
    node_betweenness <- betweenness(interaction_network, directed = FALSE, weights = E(interaction_network)$weight)
    
    # Add betweenness centrality to node data
    node_data <- data.frame(
      name = V(interaction_network)$name,
      X = mif_data$X[match(V(interaction_network)$name, mif_data$Cell_ID)],
      Y = mif_data$Y[match(V(interaction_network)$name, mif_data$Cell_ID)],
      cell_type = mif_data$cell_type[match(V(interaction_network)$name, mif_data$Cell_ID)],
      betweenness = node_betweenness
    )
    
    # Extract edge data for visualization
    network_edges <- as.data.frame(get.data.frame(interaction_network, what = "edges"))
    
    # Add spatial coordinates for edges
    network_edges <- network_edges %>%
      mutate(
        x_start = mif_data$X[match(from, mif_data$Cell_ID)],
        y_start = mif_data$Y[match(from, mif_data$Cell_ID)],
        x_end = mif_data$X[match(to, mif_data$Cell_ID)],
        y_end = mif_data$Y[match(to, mif_data$Cell_ID)]
      )
    
    # Visualize the network with betweenness centrality
    p = ggplot() +
      # Plot edges
      geom_segment(data = network_edges,
                   aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                   color = "gray", alpha = 0.5, size = 0.5) +
      # Plot nodes, sized and colored by betweenness centrality
      geom_point(data = node_data, aes(X, Y, size = betweenness, color = cell_type), alpha = 0.8) +
      scale_size_continuous(name = "Betweenness Centrality") +
      labs(title = "Interaction Network with Betweenness Centrality",
           x = "X Coordinate", y = "Y Coordinate") +
      theme_minimal()
    
    ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}
