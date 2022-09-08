chemtaxbrokewest <- function(idx = c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1),
                             type = "sd") {
################################################################################
#                                                                              # 
#            ---- Read phytoplankton data for chemtaxbrokewest ----            #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# The chemtaxbrokewest.csv or b is read that contains the pigment data for all
# measured pigments (rows = sample, col = pigment). The
# brokewest_pigment_ratios.csv or b is read that contains the initial pigments
# ratio matrix (row = taxa ratio, col 1 = species name, others = pigment ratio).
# The pigments are selected by an index of 1/0 representing the cols in the
# pigment ratio matrix that is used. This is pre-selected for testing. A
# comparison between pigments names from pigment ratio matrix and sample matrix
# is done to select the same pigments (spelling counts) and check that none is
# missing from the sample matrix. Two matrices are created for initial standard
# deviation using the formula(s):
#
# df_sd      = df_sd * 0.01 + 0.0003;
# pig_r_sd   = pig_r_init*0.1; and
# last col   = 0.05 
# 
# -or- using normalizing
# 
# df_norm    = df / row sum df
# pig_r_nrom = pig_r / row sum pig_r
# 
# for sample matrix and pigment ratio matrix, respectively. The pigment names  
# and taxa name are exported as their own variable as well.
#
# ---- INPUTS: -----------
# idx        = index of selected and deselected pigment from pigment ratio matrix
# norm       = if want to create sd matrix or normalized matrix; 
#              options: sd   = standard deviations
#                       norm = normalized to row sums
#
# ---- OUTPUTS: ----------
# df         = Matrix of samples by pigment readings
# pig_r_init = Initial pigment ratio matrix by taxa
# taxa       = Cell array of taxa names
# pigm_sel   = Cell array of selected pigment names
#
# Either:
# df_sd      = Standard deviations for df_sd (i.e. df_sd*0.01+0.0003)  
# pig_r_sd   = Standard deviations for f (i.e pig_r_init*0.1 & last col = 0.05)
#
# -or- 
#
# df_norm    = Normalized df to row sums
# pig_r_norm = Normalized pigment ratio matrix to row sums
#  
# ---- NOTES: ------------
# Original: 2010-03-21  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Mon Jun 13 18:21:17 2022)
  
  # ---- set directory ----
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/" 
  source(paste0(root,"/scripts/permcalc.R"))
  
  # ---- read data ----
  # set file paths
  # TODO: should make this a function?
  sample_filepath      <- switch(
    type,
    norm = paste0(root, raw, "CHEMTAXBROKEWestx.csv"),
    sd   = paste0(root, raw, "CHEMTAXBROKEWests.csv")
  )
  
  pig_ratio_filepath   <- switch(
    type,
    norm = paste0(root, "/scripts/brokewest_b_pigment_ratios.csv"),
    sd   = paste0(root, "/scripts/brokewest_pigment_ratios.csv"),
  )
  
  # pigment concentration in samples
  df                   <- as.matrix(read.csv(sample_filepath))
  
  # pigment ratio matrix
  temp_ratio           <- read.csv(pig_ratio_filepath)
  
  # extract pigment names
  df_pig               <- colnames(df)
  pigm_temp <- colnames(temp_ratio)[-1]
  
  colnames(df)         <- NULL # TODO: decide if it matters or not
  colnames(temp_ratio) <- NULL # TODO: decide if it matters or not
  
  # extract taxa names
  taxa                 <- temp_ratio[,1]
  
  # ---- TODO: remove later ----
  # idx <- c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)
  
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
  df  <- df[, df_pig_idx]
  
  if (type != "norm") {
    # ---- set standard deviation values ----
    # TODO: make this part of function
    df_sd                      <- df * 0.01 + 0.0003
    pig_r_sd                   <- pig_r_init * 0.1
    pig_r_sd[, ncol(pig_r_sd)] <- 0.005 # set chlor-a sd to 0.005
    
    # ---- return list of variables created ----
    result <-
      list(
        df         = df,
        df_sd      = df_sd,
        pig_r_init = pig_r_init,
        pig_r_sd   = pig_r_sd ,
        taxa       = taxa,
        pigm_sel   = pigm_sel
      )
  } else {
    # ---- normalize df and pigment ratio to row sums ----
    df_row_sum    <- as.matrix(apply(df, 1, sum, na.rm = TRUE))
    pig_r_row_sum <-
      as.matrix(apply(pig_r_init, 1, sum, na.rm = TRUE))
    
    df_norm       <- df / pracma::repmat(df_row_sum, 1, ncol(df))
    pig_r_norm    <-
      pig_r_init / pracma::repmat(pig_r_row_sum, 1, ncol(pig_r_init))
    
    # ---- return list of variables created ----
    result <-
      list(
        df         = df,
        df_norm    = df_norm,
        pig_r_init = pig_r_init,
        pig_r_norm = pig_r_norm ,
        taxa       = taxa,
        pigm_sel   = pigm_sel
      )
  }
}



