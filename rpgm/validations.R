

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

  required_properties <- c("sample_id", #"tissue_class_label",
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

validate_bivar = function(config){
  has_bivar = any(c("kcross", "gcross", "kcross_exactCSR", "gcross_perm") %in% config$metrics)
  if(has_bivar){
    if(inherits(config$variables$bivar_pairs, "list")){
      out = lapply(config$variables$bivar_pairs, function(x){
        #if someone wants marker v all:
        if("all" %in% x){
          w_all = which(x == "all")
          #okay check which one is all and expand
          if(length(w_all) == 2){
            pairs = expand.grid(marker1 = config$variables$markers,
                                marker2 = config$variables$markers)
          } else {
            if(w_all == 1){
              pairs = expand.grid(marker1 = config$variables$markers,
                                  marker2 = x[2])
            } else {
              pairs = expand.grid(marker1 = x[1],
                                  marker2 = config$variables$markers)
            }
          }
          
        } else {
          pairs = data.frame(marker1 = x[1],
                             marker2 = x[2])
        }
        return(pairs)
      }) %>%
        do.call(dplyr::bind_rows, .) %>%
        dplyr::filter(marker1 != marker2)
      config$variables$bivar_pairs = lapply(seq(nrow(out)), function(x) as.character(out[x,]))
    } else {
      stop("Error: Bivariate pairs must be format - ['marker', 'marker']")
    }
  }
  return(config)
}
