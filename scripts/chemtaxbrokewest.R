chemtaxbrokewest <- function() {
################################################################################
#                                                                              # 
#            ---- Read phytoplankton data for chemtaxbrokewest ----            #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# The chemtaxbrokewest.csv is read that contains the pigment data for all
# measured pigments (rows = sample, col = pigment). The 
# brokewest_pigment_ratios.csv is read that contains the initial pigments ratio
# matrix (row = taxa ratio, col 1 = species name, others = pigment ratio). The
# pigments are selected by an index of 1/0 representing the cols in the pigment 
# ratio matrix that is used. This is pre-selected for testing. A comparison 
# between pigments names from pigment ratio matrix and sample matrix is done to
# select the same pigments (spelling counts) and check that none is missing from
# the sample matrix. Two matrices are created for initial standard deviation 
# using the formula(s):
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
  
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/"
  
  # ---- CHEMTAXBROKEWests  ----
  # pigment concentration in samples
  # TODO: should make this a function?
  sample_filepath     <- paste0(root, raw, "CHEMTAXBROKEWests.csv")
  df           <- as.matrix(read.csv(sample_filepath))
  
  # extract pigment names
  df_pig              <- colnames(df)
  colnames(df) <- NULL
 
  # ---- index for selected pigments of pigment ratios matrix ----
  # TODO: make input for this 
  idx <- c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)
  
  # ---- read pigment ratios ----
  pig_ratio_filepath       <- paste0(root,"/scripts/brokewest_pigment_ratios.csv")
  temp_ratio               <- read.csv(pig_ratio_filepath)
  
  # initialize ratio matrix by removing taxa and pigment names
  pig_r_init           <- as.matrix(temp_ratio[,-1])[, which(idx==1)] 
  colnames(pig_r_init) <- NULL # TODO: decide if it matters or not
  
  # extract pigment names
  pigm_sel                 <- colnames(temp_ratio)[-1][which(idx==1)] 
  
  # extract taxa names
  taxa                     <- temp_ratio[,1]

  # ---- find columns that match and rearrange ----
  source(paste0(root,"/scripts/permcalc.R"))
  df_pig_idx <- permcalc(pigm_sel, df_pig)
  
  # filter columns for pigment that match selected pigments
  df  <- df[, df_pig_idx]
  
  # ---- set standard deviation values ----
  df_sd                       <- df * 0.01 + 0.0003
  pig_r_sd                   <- pig_r_init * 0.1
  pig_r_sd [, ncol(pig_r_sd )] <- 0.005 # set chlor-a sd to 0.005

  # return list of all variables created
  result <-
    list(
      df      = df,
      df_sd   = df_sd,
      pig_r_init = pig_r_init,
      pig_r_sd    = pig_r_sd ,
      taxa           = taxa,
      pigm_sel       = pigm_sel
    )
  
}



