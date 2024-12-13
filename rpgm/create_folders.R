createfolders <- function(yaml){ 
  output_path <- yaml$paths$output
  if(!dir.exists(output_path)){
    dir.create(output_path)
  }
  if(!dir.exists(file.path(output_path, "figures"))){
    dir.create(file.path(output_path, "figures"))
  }
  if(!dir.exists(file.path(output_path, "metrics"))){
    dir.create(file.path(output_path, "metrics"))
  }
  
  figures_path <- file.path(output_path, "figures")
  if(!dir.exists(file.path(figures_path, "barplot"))){
    dir.create(file.path(figures_path, "barplot"))
  }
  if(!dir.exists(file.path(figures_path, "metrics"))){
    dir.create(file.path(figures_path, "metrics"))
  }
  if(!dir.exists(file.path(figures_path, "point"))){
    dir.create(file.path(figures_path, "point"))
  }
  
  #prep metrics folders
  tmp = lapply(yaml$metrics, function(m){
    if(!dir.exists(file.path(figures_path, "metrics", m))){
      dir.create(file.path(figures_path, "metrics", m))
    }
    if(!dir.exists(file.path(output_path, "metrics", m))){
      dir.create(file.path(output_path, "metrics", m))
    }
  })
  
}
