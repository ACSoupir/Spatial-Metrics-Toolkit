#packages
renv::activate()

install.packages(
  c("spatstat", "dbscan", "dixon", #analyses
    "magrittr", "data.table", "dplyr", "tibble", #data processing
    "spatialTIME", #example data
    "optparse" #allow Rscript to run as command line tool
    ),
  Ncpus = 4
)

library(dplyr)
library(data.table)
library(spatialTIME)
#prepping data from the spatialTIME package
spatial_list = example_spatial
#get columns we want
sample_id = "deidentified_sample"
markers = c("CD3..FOXP3.", "CD3..CD8.", "CD3..CD8..FOXP3.",
            "CD3..PD1.", "CD3..PD.L1.", "CD8..PD1.", "CD3..CD8..PD.L1.",
            "CD3..Opal.570..Positive", "CD8..Opal.520..Positive",
            "PD1..Opal.650..Positive", "FOXP3..Opal.620..Positive",
            "PDL1..Opal.540..Positive")
tissue_col = "Classifier.Label"
#iterate over each sample
spatial_list2 = lapply(spatial_list, function(spat){
  df = spat %>%
    mutate(x = (XMin + XMax)/2, #x and y are ranges
           y = (YMin + YMax)/2) %>%
    select(any_of(c(sample_id, "x", "y", markers, tissue_col))) #select our columns
  num_pos = sample(15:150, size = length(markers))
  names(num_pos) = markers
  for(marker in markers){
    cells = sample(1:nrow(df), size = num_pos[marker])
    df[cells,marker] = 1
  }
  return(df)
})

#export to data/per-cell/
tmp = lapply(spatial_list2, function(spat){
  path = paste0("data/per-cell/",
                unique(spat[[sample_id]]),
                ".csv.gz")
  fwrite(spat, file = path, compress = 'gzip')
})
