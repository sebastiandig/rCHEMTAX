load_data <- function(.file = NULL, .pig_r = NULL, .pig_sel = NULL,
                      .file_sd = NULL, .pig_r_sd = NULL, .clust_col = NULL,
                      verbose = TRUE) {
  
  
  # TODO: make function to read in data set, and pigment ratio matrix,
  #       select pigments, and if want to use generic SD or have specific info
  #       for both pigment ratio matrix, data set
  
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/" 
  .file      = paste0(root, raw, "CHEMTAXBROKEWestx.csv")

  .pig_file     = paste0(root, "/scripts/brokewest_b_pigment_ratios.csv")
  .pig_sel   = NULL
  .file_sd   = NULL
  .pig_r_sd  = NULL
  .clust_col = NULL
  verbose    = TRUE
  
  
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
    csv  = read.csv(.file),
    xlsx = readxl::read_xlsx(.file),
    xls  = readxl::read_xls(.file),
    "NA" = stop("file extension is not .csv, .xlsx, or .xls")
  )
  
  temp_ratio <- switch(str_extract(.pig_file, "csv|xlsx|xls"),
                csv  = read.csv(.pig_file),
                xlsx = readxl::read_xlsx(.pig_file),
                xls  = readxl::read_xls(.pig_file),
                "NA" = stop("file extension is not .csv, .xlsx, or .xls")
  )
   
  # clusters <-
  #   switch(.clust_col,
  #  !is.null(.clust_col) = df1[[.clust_col]] 
  # )
  
  # .clust_col can either be number of name
  
  clusters <- NULL

  clusters <- switch(!is.null(.clust_col),
    "TRUE" = df1[[.clust_col]],
    "FALSE" = matrix(1, nrow = nrow(df1), ncol = 1)
  ) 
  
  if (!is.null(.clust_col)) {
    clusters <- df1[[.clust_col]]
  } 

  pigm_temp <- colnames(temp_ratio)[-1]
  
  # extract taxa names
  stopifnot("First column in pigment ratio matrix needs to be a string of taxa!" 
            = is.character(temp_ratio[,1]))
  taxa                 <- temp_ratio[,1]
}
