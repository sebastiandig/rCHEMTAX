chemtaxbrokewest_b <- function(idx = c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)) {
################################################################################
#                                                                              # 
#           ---- Read phytoplankton data for chemtaxbrokewestb ----            #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# The CHEMTAXBROKEWestx.csv is read that contains the pigment data for all
# measured pigments (rows = sample, col = pigment). The
# brokewest_b_pigment_ratios.csv is read that contains the initial pigments
# ratio matrix (row = taxa ratio, col 1 = species name, others = pigment ratio).
# The pigments are selected by an index of 1/0 representing the cols in the
# pigment ratio matrix that is used. This is pre-selected for testing. A
# comparison between pigments names from pigment ratio matrix and sample matrix
# is done to select the same pigments (spelling counts) and check that none is
# missing from the sample matrix. Two matrices are created for initial standard
# deviation using the formula(s):
# 
# df_norm    = df / row sum df
# pig_r_nrom = pig_r / row sum pig_r
# 
# ---- INPUTS: -----------
# idx        = index of selected and deselected pigment from pigment ratio matrix
#
# ---- OUTPUTS: ----------
# df         = Matrix of samples by pigment readings
# df_norm    = Normalized df to row sums
# pig_r_init = Initial pigment ratio matrix by taxa
# pig_r_norm = Normalized pigment ratio matrix to row sums
# taxa       = Cell array of taxa names
# pigm_sel   = Cell array of selected pigment names 
#
# ---- NOTES: ------------
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (2022-09-08 18:59:47)  
  # ---- load library ----
  library("pracma")
  
  # ---- set directory ----
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/"
  source(paste0(root, "/scripts/permcalc.R"))
  # TODO: matfactuvw
  source(paste0(root, "/scripts/matfactuvw.R"))
  
  # ---- read data ----
  # CHEMTAXBROKEWestx pigment concentration
  # TODO: should make this a function?
  sample_filepath <- paste0(root, raw, "CHEMTAXBROKEWestx.csv")
  df              <- as.matrix(read.csv(sample_filepath))
  
  # pigment ratio matrix
  pig_ratio_filepath   <-
    paste0(root, "/scripts/brokewest_b_pigment_ratios.csv")
  temp_ratio           <- read.csv(pig_ratio_filepath)
  
  # extract pigment names
  df_pig               <- colnames(df)
  pigm_temp <- colnames(temp_ratio)[-1]
  
  colnames(df)         <- NULL # TODO: decide if it matters or not
  colnames(temp_ratio) <- NULL # TODO: decide if it matters or not
  
  # extract taxa names
  taxa                 <- temp_ratio[, 1]
  
  # ---- index for selected pigments of pigment ratios matrix ----
  if (is.null(idx) | length(idx) != length(pigm_temp)) {
    message(
      paste0(
        "Either an index was not supplied, was 0, or exceeded the number of pigments!\n",
        "\nHere are a list of pigment names:\n"
      )
    )
    message(paste(pigm_temp, "\n"))
    message(paste0("Give a 1 or a 0 for each name\n",
                   "(1 = selected; 0 = not selected)\n"))
    idx <- readline(prompt = paste0("Should be a length of ",
                                    length(pigm_temp),
                                    " (ex 1 0 1 1 0):\n "))
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
  pig_r_init           <-
    as.matrix(temp_ratio[, -1])[, which(idx == 1)]
  
  # extract selected pigment names
  pigm_sel             <- pigm_temp[which(idx == 1)]
  
  # ---- find columns that match and rearrange ----
  df_pig_idx <- permcalc(pigm_sel, df_pig)
  
  # filter columns for pigment that match selected pigments
  df  <- df[, df_pig_idx]
  
  # ---- normalize df and pigment ratio to row sums ----
  df_row_sum    <- as.matrix(apply(df, 1, sum, na.rm = TRUE))
  pig_r_row_sum <-
    as.matrix(apply(pig_r_init, 1, sum, na.rm = TRUE))
  
  df_norm       <- df / pracma::repmat(df_row_sum, 1, ncol(df))
  pig_r_norm    <- pig_r_init / pracma::repmat(pig_r_row_sum, 1, ncol(b1))
  
  # ---- return list of variables created ----
  result <-
    list(
      df         = df,
      df_norm = df_norm,
      pig_r_init = pig_r_init,
      pig_r_norm = pig_r_norm ,
      taxa       = taxa,
      pigm_sel   = pigm_sel
    )
  
}