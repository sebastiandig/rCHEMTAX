amatfactsd <- function(.df, .df_sd, .pig_r){
################################################################################
#                                                                              # 
#               Calculate Factor a in min(x-ab) using least square             #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# Solves for a
# brokewest -> chemtaxbrokewest() -> line 57 -> nnmatfactsd -> amatfactsd
# 
# ---- INPUTS: -----------
# df           = Matrix to be factored as taxa_amt_min*pig_r
# sdx          = Standard deviations of df
# pig_r        = Value of pigment ratios
#
# ---- OUTPUTS: ----------
# taxa_amt_min = Result left factor of df
#
# ---- NOTES: -------------
# Original: 2010-03-14  Matalb7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Thu Jun 02 18:30:06 2022)
library("pracma")
  
  # ---- row numbers of the data set and pigment ratios matrix ----
  ns <- dim(.df)[1]
  nt <- dim(.pig_r)[1]
  
  # ---- initialize empty matrix ----
  taxa_amt_min <- matrix(0, ns, nt)
  
  # ---- loop through rows of df ----
  for (j in seq(ns)) {
    
    # ---- setup row values ----
    sdx_j   <- pracma::repmat(as.matrix(.df_sd[j,]), 1, nt)
    df_j    <- as.vector(t((.df[j,] / .df_sd[j,])))
    pig_r_j <- t(.pig_r) / sdx_j
    
    # ---- non-negative least linear square minimization ----
    temp    <- pracma::lsqnonneg(pig_r_j, df_j)
    
    # ---- results transposed ----
    taxa_amt_min[j,] <- t(temp$x)
  }
  
  taxa_amt_min
  
}
          