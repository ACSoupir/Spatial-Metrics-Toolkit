

validate_property <- function(data, property) {
  if (is.null(data[[property]])) {
    stop(paste("Error: The required property", property, 
               "is missing in the config file."))
  }
}

validate_property_is_a_list <- function(data, property) {
  if (!((is.list(data[[property]]) || is.vector(data[[property]])) 
        && any(sapply(data[[property]], is.character)))) {
    stop(paste("Error: The property", property, 
               "does not contain a list of strings."))
  }
}

#Validate config data
validate_config <- function(config) {

  validate_property(config, "variables")

  required_properties <- c("sample_id", "tissue_class_label",
                           "markers", "radii_range")
  for (property in required_properties) {
    validate_property(config$variables, property)
  }

  validate_property(config, "metrics")
  validate_property_is_a_list(config, "metrics")

  validate_property(config, "paths")
  required_properties <- c("spatial", "output", "sample")
  for (property in required_properties) {
    validate_property(config$paths, property)
  }

}

validate_df_columns <- function(df, columns) {
  if (!all(columns %in% colnames(df))) {
    missing_columns <- columns[!columns %in% colnames(df)]
    stop(paste("Error: The following required columns are missing:", paste(missing_columns, collapse = ", ")))
  }

}
