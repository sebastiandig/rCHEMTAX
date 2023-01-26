bmatfactsd <- function(.df, .df_sd, .taxa_amt, .pig_r, .pig_r_sd) {
  ################################################################################
  #                                                                              # 
  #                  Calculate factor b in min(x-ab, b-b0) using                 #
  #                       Non-negative Least Linear Squares                      #
  #                                                                              #
  ################################################################################
  # ---- DESCRIPTION: ------
  # Solve for b
  # This calculates the pigment ratio matrix using the taxa contribution matrix 
  # from the results of `nnmatfactsd` algorithm. The pigment ratio matrix is `b` 
  # from min||x-ab, b-b0||.
  #
  # ---- INPUTS: -----------
  # df        = Matrix to be factored as a*b
  # df_sd     = Standard deviations of df
  # taxa_amt  = result of amatfactsd(x,sdx,b,sdb)
  # pig_r     = Initial estimate for b
  # pig_r_sd  = Standard deviations of pig_r?
  #        
  # ---- OUTPUTS: ----------
  # pig_r_new = Result left factor of x
  # 
  # ---- NOTES: ------------
  # Original: 2010-03-14  Matalb7  W.Whiten
  #
  # ---- REFERENCES(s): ----
  #
  # ---- AUTHOR(s): --------
  # Sebastian Di Geronimo (Thu Jun 02 19:54:42 2022)
  
  # ---- remove after testing ----
  # .df       = s
  # .df_sd    = ssd
  # .taxa_amt = taxa_amt
  # .pig_r    = f
  # .pig_r_sd = fsd

  
  # ---- load library ----
  library("pracma")
  
  # ---- dimensions of pigment ratio matrix ----
  .pig_r_row <- dim(.pig_r)[1] # row
  .pig_r_col <- dim(.pig_r)[2] # col
  
  # index where b0 does not equal 0
  indx      <- .pig_r != 0
  
  # create empty matrix of 0s with nt and np dimenstions
  pig_r_new <- matrix(0, .pig_r_row, .pig_r_col)
  
  # ---- iterate through each row in pigment ratio then calc b ----
  for (j in seq(.pig_r_col)) {
    idx <- indx[,j]
    
    # ---- setup variables for lsqnonneg  ----
    # take columns in taxa_j where has the pigment in df_sd[,j], divide by sd
    taxa_by_df_sd <- .taxa_amt[,idx] / pracma::repmat(as.matrix(.df_sd[,j]),1, sum(idx))
    # pig_r_sd_diag <- diag(1 / .pig_r_sd[idx,j])
    
    if (length(.pig_r_sd[idx,j]) == 1) {
      pig_r_sd_diag <- diag(1 / .pig_r_sd[idx,j], nrow = length(1))
    } else {
      pig_r_sd_diag <- diag(1 / .pig_r_sd[idx,j])
    }
    
    # vars pigment column j and divide by its sd
    df_by_sd      <- .df[,j] / .df_sd[,j]
    pig_by_sd     <- .pig_r[idx,j] / .pig_r_sd[idx,j]
    
    # combine data into matrix and vector for 
    taxa_norm     <- as.matrix(rbind(taxa_by_df_sd, pig_r_sd_diag))
    df_pig_norm   <- as.vector(c(df_by_sd,pig_by_sd))
    
    # ---- lsqnonneg (formula min||x - a*b||) ----
    # gives `a` with taxa_norm
    b_j <- pracma::lsqnonneg(
                              taxa_norm,
                              df_pig_norm
                             )

    pig_r_new[idx,j] <- t(b_j$x) 
  }
  
  # ---- return values ----
  pig_r_new
}