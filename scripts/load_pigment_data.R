load_data <- function(.file = NULL, .pig_r = NULL, .pig_sel = NULL,
                      .file_sd = NULL, .pig_r_sd = NULL, .clust_col = NULL,
                      verbose = TRUE) {
  
  
  # TODO: make function to read in data set, and pigment ratio matrix,
  #       select pigments, and if want to use generic SD or have specific info
  #       for both pigment ratio matrix, data set
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/" 
  .file = paste0(root, raw, "CHEMTAXBROKEWestx.csv")

  .pig_r = NULL
  .pig_sel = NULL
  .file_sd = NULL
  .pig_r_sd = NULL
  .clust_col = NULL
  verbose = TRUE
  
  
  library("ggplot2")
  library("tibble")
  library("tidyr")
  library("readr")
  library("purrr")
  library("dplyr")
  library("stringr")
  library("forcats")
  library("lubridate")
  library("glue")
  library("fs")
  library("magrittr")
  # library("broom") # optional
  
  assertthat::is.string(.file)
  

  
  ?setClass
  ??assertthat
  ?methods
  
  df1 <- switch(str_extract(.file, "csv|xlsx|xls"),
    csv = read.csv(.file),
    xlsx = readxl::read_xlsx(.file),
    xls = readxl::read_xls(.file)
  )
  
  # clusters <-
    switch(.clust_col,
   !is.null(.clust_col) = df1[[.clust_col]] 
  )
  
  clusters <- NULL
  if (!is.null(.clust_col)) {
    clusters <- df1[[.clust_col]]
  } 



  
  
  
}