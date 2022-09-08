chemtaxbrokewest <- function(idx = c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)) {
################################################################################
#                                                                              # 
#            ---- Read phytoplankton data for chemtaxbrokewest ----            #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# The chemtaxbrokewest.csv is read that contains the pigment data for all
# measured pigments (rows = sample, col = pigment). The
# brokewest_pigment_ratios.csv is read that contains the initial pigments
# ratiomatrix (row = taxa ratio, col 1 = species name, others = pigment ratio).
# The pigments are selected by an index of 1/0 representing the cols in the
# pigment ratio matrix that is used. This is pre-selected for testing. A
# comparison between pigments names from pigment ratio matrix and sample matrix
# is done to select the same pigments (spelling counts) and check that none is
# missing from the sample matrix. Two matrices are created for initial standard
# deviation using the formula(s):
#
# df_sd = df_sd * 0.01 + 0.0003;
# pig_r_sd  = pig_r_init*0.1; and
# last col     = 0.05 
# 
# for sample matrix and pigment ratio matrix, respectively. The pigment names  
# and taxa name are exported as their own variable as well.
#
# ---- INPUTS: -----------
# NA
#
# ---- OUTPUTS: ----------
# df   = Matrix of samples by pigment readings
# df_sd   = Standard deviations for df_sd (i.e. df_sd*0.01+0.0003)
# pig_r_init = Initial pigment ratio matrix by taxa
# pig_r_sd     = Standard deviations for f (i.e pig_r_init*0.1 & last col = 0.05)
# taxa           = Cell array of taxa names
# pigm_sel       = Cell array of selected pigment names
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
  
  
  # ---- read data ----
  # pigment concentration in samples
  # TODO: should make this a function?
  sample_filepath <- paste0(root, raw, "CHEMTAXBROKEWests.csv")
  df              <- as.matrix(read.csv(sample_filepath))
  
  # pigment ratio matrix
  pig_ratio_filepath   <- paste0(root, "/scripts/brokewest_pigment_ratios.csv")
  temp_ratio           <- read.csv(pig_ratio_filepath)
  
  # extract pigment names
  df_pig       <- colnames(df)
  colnames(df) <- NULL
 
  # extract pigment names
  pigm_temp <- colnames(temp_ratio)[-1]
  
  
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
  colnames(pig_r_init) <- NULL # TODO: decide if it matters or not
  
  # extract pigment names
  pigm_sel             <- colnames(temp_ratio)[-1][which(idx == 1)] 
  
  # extract taxa names
  taxa                 <- temp_ratio[,1]

  # ---- find columns that match and rearrange ----
  source(paste0(root,"/scripts/permcalc.R"))
  df_pig_idx <- permcalc(pigm_sel, df_pig)
  
  # filter columns for pigment that match selected pigments
  df  <- df[, df_pig_idx]
  
  # ---- set standard deviation values ----
  # TODO: make this part of function
  df_sd                       <- df * 0.01 + 0.0003
  pig_r_sd                    <- pig_r_init * 0.1
  pig_r_sd[, ncol(pig_r_sd )] <- 0.005 # set chlor-a sd to 0.005

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
  
}



