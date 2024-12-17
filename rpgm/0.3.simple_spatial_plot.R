
#' Title: simple_spatial_plot
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param x x axis
#' @param y y axis
#' @param pos_is_cell_positive expect 0 or 1
#'
#' @return ggplot2 plot object
#' @examples
#' # Example usage
#' simple_spatial_plot(data_base, "CD3..FOXP3.", "Classifier.Label")
#'
#' @export
simple_spatial_plot <- function(data_base, marker, shape) {
  #data_base <- read.csv(filepath, header = TRUE, sep = ",")
  #p1 = simple_spatial_plot(data_base, colname)
  
  plot_data = cbind(data_base$x, data_base$y, data_base[,marker], data_base[,shape])
  
  colnames(plot_data) = c("x", "y", marker, shape)
  plot_data = as.data.frame(plot_data)
  
  plot_data[1:2] <- lapply(plot_data[1:2], as.numeric) # forcing the x and y as numeric
  
  plot_data[which(plot_data[,marker] == 1), marker] = "present"
  plot_data[which(plot_data[,marker] == 0), marker] = "absent"
  p = ggplot(plot_data, aes(x = x, y = y, color = !!sym(marker), shape = !!sym(shape))) + geom_point()
  return(p)
}


#' Title: simple_spatial_plot
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param x x axis
#' @param y y axis
#' @param pos_is_cell_positive expect 0 or 1
#'
#' @return ggplot2 plot object
#' @examples
#' # Example usage
#' simple_spatial_plot(data_base, "CD3..FOXP3.", "Classifier.Label")
#'
#' @export
simple_spatial_plot <- function(data_base, marker, shape, merge_marker_shape) {
  #data_base <- read.csv(filepath, header = TRUE, sep = ",")
  #p1 = simple_spatial_plot(data_base, colname)
  
  data_base[which(data_base[,marker] == 1), marker] = "pos"
  data_base[which(data_base[,marker] == 0), marker] = "neg"
  if (merge_marker_shape) {
    plot_data = cbind(data_base$x, data_base$y, paste0(data_base[,marker], "_", data_base[,shape]))
    maker_shape =  paste0(marker, "_", shape)
    colnames(plot_data) = c("x", "y", maker_shape)
    plot_data = as.data.frame(plot_data)
    
    plot_data[1:2] <- lapply(plot_data[1:2], as.numeric) # forcing the x and y as numeric
    
    p = ggplot(plot_data, aes(x = x, y = y, color = !!sym(maker_shape))) + geom_point()
    return(p)
    
  } else {
    plot_data = cbind(data_base$x, data_base$y, data_base[,marker], data_base[,shape])
    
    colnames(plot_data) = c("x", "y", marker, shape)
    plot_data = as.data.frame(plot_data)
    
    plot_data[1:2] <- lapply(plot_data[1:2], as.numeric) # forcing the x and y as numeric
    
    p = ggplot(plot_data, aes(x = x, y = y, color = !!sym(marker), shape = !!sym(shape))) + geom_point()
    return(p)
    
  }
}

#' Title: Generate a bunch of plots
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#'
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generatePlotAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generatePlotAsPDF <- function(yamlfile, csvfilePath) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/", clean_name, "_", marker, ".pdf")
    p = simple_spatial_plot(data_base, marker, categories, FALSE)
    ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}


#' Title: Generate a bunch of plots
#'
#' Takes in a csv file expecting x and y with a value
#'
#' @param yamlfile yaml file containing the marker list and cell type column
#' @param csvfilepath of the matrix
#'
#' @return output folder with pdf
#' @examples
#' # Example usage
#' generatePlotWithCelltypeMergeAsPDF("config.yml", "data/mIF_per-cell/TMA3_[9,K].tif.csv.gz")
#'
#' @export
generatePlotWithCelltypeMergeAsPDF <- function(yamlfile, csvfilePath) {
  
  clean_name <- gsub("\\.gz$", "", basename(csvfilePath))
  clean_name <- gsub("\\.csv$", "", clean_name)
  clean_name <- gsub("\\.txt$", "", clean_name)
  
  data_base <- read.csv(filepath, header = TRUE, sep = ",")
  outputFilePath = read_yaml(yamlfile)$paths$output
  categories = read_yaml(yamlfile)$variables$tissue_class_label
  
  markers = read_yaml(yamlfile)$variables$markers
  for (marker in markers) {
    outputPDFname = paste0(outputFilePath, "/", clean_name, "_", marker, ".pdf")
    p = simple_spatial_plot(data_base, marker, categories, TRUE)
    ggsave(outputPDFname, plot = p, width = 8, height = 6)
  }
}
