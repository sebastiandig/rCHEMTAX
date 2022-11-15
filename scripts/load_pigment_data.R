load_data <- function(.file = NULL, .pig_file = NULL, idx = NULL,
                      .file_sd = NULL, .pig_r_sd = NULL, .clust_col = NULL, 
                      type = "sd",
                      verbose = TRUE) {
  
  
  # TODO: make function to read in data set, and pigment ratio matrix,
  #       select pigments, and if want to use generic SD or have specific info
  #       for both pigment ratio matrix, data set
  
  root <- rprojroot::find_rstudio_root_file()
  # raw  <- "/data/raw/" 
  
  library("stringr")
  
  source(paste0(root,"/scripts/permcalc.R"))
  
  # .file      = paste0(root, raw, "CHEMTAXBROKEWests.csv")
  # .pig_file     = paste0(root, "/scripts/brokewest_pigment_ratios.csv")
  # idx = c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)
  # .file_sd   = NULL
  # .pig_r_sd  = NULL
  # .clust_col = NULL
  # verbose    = TRUE
  # type = "sd"



  # library("ggplot2")
  # library("tibble")
  # library("tidyr")
  # library("readr")
  # library("purrr")
  # library("dplyr")

  # library("forcats")
  # library("lubridate")
  # library("glue")
  # library("fs")
  # library("magrittr")
  # library("broom") # optional
  
  # assertthat::is.string(.file)
  
  # ---- check and load files ----
  stopifnot("Data filepath was not given." = !is.null(.file))
  stopifnot("Pigment ratio filepath was not given." = !is.null(.pig_file))
  
  df <- switch(
    stringr::str_extract(.file, "csv|xlsx|xls"),
    csv  = read.csv(.file),
    xlsx = readxl::read_xlsx(.file),
    xls  = readxl::read_xls(.file),
    stop("Data filepath is extension is not .csv, .xlsx, or .xls.")
  )
 
  temp_ratio <- switch(
    stringr::str_extract(.pig_file, "csv|xlsx|xls"),
    csv  = read.csv(.pig_file),
    xlsx = readxl::read_xlsx(.pig_file),
    xls  = readxl::read_xls(.pig_file),
    stop("Pigment ratio filepath extension is not .csv, .xlsx, or .xls.")
  )

  # TODO: add switches for sd of pigment and samples
  
  # ---- extract clusters info ----
  # .clust_col can either be number or name of column with clusters info
  if (!is.null(.clust_col)) {
    clusters <- subset(df, select = .clust_col)
    
    if (is.character(.clust_col))
      .clust_col <- which(colnames(df) == .clust_col)
    
    df       <- df[-.clust_col]
    
  } else {
    
    clusters <- matrix(1, nrow = nrow(df), ncol = 1)
    
  }
  
  # ---- extract names ----
  # taxa names
  stopifnot("First column in pigment ratio matrix needs to be a string of taxa!" 
            = is.character(temp_ratio[,1]))
  taxa                 <- temp_ratio[,1]
  
  # pigment names
  df_pig               <- colnames(df)
  pigm_temp            <- colnames(temp_ratio)[-1]
  
  colnames(df)         <- NULL # TODO: decide if it matters or not
  colnames(temp_ratio) <- NULL # TODO: decide if it matters or not
  
  # ---- index for selected pigments of pigment ratios matrix ----
  if (is.null(idx) | length(idx) != length(pigm_temp)) {
    message(paste0("Either an index was not supplied, was 0, or exceeded the number of pigments!\n",
                   "\nHere are a list of pigment names:\n"))
    message(paste(pigm_temp, "\n"))
    message(paste0("Give a 1 or a 0 for each name\n",
                   "(1 = selected; 0 = not selected)\n"))
    idx <- readline(prompt = paste0(
      "Should be a length of ", 
      length(pigm_temp), 
      " (ex 1 0 1 1 0):\n ")
    )
    idx <- as.integer(strsplit(idx, " ")[[1]])
    if (length(idx) != length(pigm_temp)) {
      idx <-  NULL
      stop("Input does not match length of pigments") 
    }
    message(paste0("\nVerify selected pigments are correct:\n"))
    cat(pigm_temp[which(idx == 1)])
  }
  
  # ---- read pigment ratios ----
  # initialize ratio matrix by removing taxa and pigment names
  pig_r_init           <- as.matrix(temp_ratio[,-1])[, which(idx == 1)] 
  
  # extract selected pigment names
  pigm_sel             <- pigm_temp[which(idx == 1)]
  
  # ---- find columns that match and rearrange ----
  
  df_pig_idx <- permcalc(pigm_sel, df_pig)
  
  # filter columns for pigment that match selected pigments
  df  <- as.matrix(df[, df_pig_idx])
  
  # ---- set variations in df and pig_r ----
  # df
  if (is.null(.file_sd) & type != "norm") {
    # ---- standard deviation ----
    df_sd         <- df * 0.01 + 0.0003
    
  } else if (is.null(.file_sd)) {
    # ---- normalize df and pigment ratio to row sums ----
    df_row_sum    <- as.matrix(apply(df, 1, sum, na.rm = TRUE))
    df_sd         <- df / pracma::repmat(df_row_sum, 1, ncol(df))
    
  } else {
    df_sd         <-  .file_sd
  }
  
  # pig_r
  if (is.null(.pig_r_sd) & type != "norm") {
    # ---- standard deviation ----
    pig_r_sd                   <- pig_r_init * 0.1
    pig_r_sd[, ncol(pig_r_sd)] <- 0.005 # set chlor-a sd to 0.005
    
  } else if (is.null(.pig_r_sd)) {
    # ---- normalize df and pigment ratio to row sums ----
    pig_r_row_sum              <-
      as.matrix(apply(pig_r_init, 1, sum, na.rm = TRUE))
    pig_r_sd                   <-
      pig_r_init / pracma::repmat(pig_r_row_sum, 1, ncol(pig_r_init))
    
  } else {
    pig_r_sd                   <-  .pig_r_sd
  }

  # ---- return list of variables created ----
  result <-
    list(
      df         = df,
      df_sd      = df_sd,
      pig_r_init = pig_r_init,
      pig_r_sd   = pig_r_sd ,
      taxa       = taxa,
      pigm_sel   = pigm_sel,
      clusters   = clusters
    )
  
  # # ---- set standard deviation values ----
  # if (type != "norm") {
  #   # TODO: make this part of function
  #   df_sd                      <- df * 0.01 + 0.0003
  #   pig_r_sd                   <- pig_r_init * 0.1
  #   pig_r_sd[, ncol(pig_r_sd)] <- 0.005 # set chlor-a sd to 0.005
  #   
  #   # ---- return list of variables created ----
  #   result <-
  #     list(
  #       df         = df,
  #       df_sd      = df_sd,
  #       pig_r_init = pig_r_init,
  #       pig_r_sd   = pig_r_sd ,
  #       taxa       = taxa,
  #       pigm_sel   = pigm_sel,
  #       clusters   = clusters
  #     )
  # } else {
  #   # ---- normalize df and pigment ratio to row sums ----
  #   df_row_sum    <- as.matrix(apply(df, 1, sum, na.rm = TRUE))
  #   pig_r_row_sum <-
  #     as.matrix(apply(pig_r_init, 1, sum, na.rm = TRUE))
  #   
  #   df_norm       <- df / pracma::repmat(df_row_sum, 1, ncol(df))
  #   pig_r_norm    <-
  #     pig_r_init / pracma::repmat(pig_r_row_sum, 1, ncol(pig_r_init))
  #   
  #   # ---- return list of variables created ----
  #   result <-
  #     list(
  #       df         = df,
  #       df_sd      = df_sd,
  #       pig_r_init = pig_r_init,
  #       pig_r_norm = pig_r_norm ,
  #       taxa       = taxa,
  #       pigm_sel   = pigm_sel,
  #       clusters   = clusters
  #     )
  # }
  
  return(result)
  
  # ---- end ----
  }
