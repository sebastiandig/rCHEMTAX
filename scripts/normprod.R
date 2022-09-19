normprod <- function(.df,.taxa_amt,.pig_r, .col = NULL) {
################################################################################
#                                                                              # 
#                Normalize matrix product on select col of s & f               #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: -------
# After factor analysis using min||x - a * b||, results are normalized. 
#  
# The selected column needs to be the same for the data and the pigment ratio 
# matrix. If the data does not have the same number of columns, will throw and 
# error.
#
# ---- INPUTS: -----------
# df = Original product matrix
# taxa_amt = Left factor of df
# pig_r = Right factor of df
#
# ---- OUTPUTS: ----------
# ss = Matrix .df normalized to last column
# cc = Left factor after scaling
# ff = Right factor after scaling
# rms = Root mean square of ss-cc*ff
#
# ---- NOTES: ------------
# Original: 2010-03-28  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sun Jun 19 14:53:10 2022)

  # ---- load library ----
  library("pracma")
  
  # ---- extract dimension information ----
  df_row <- dim(.df)[1] # row of df
  df_col <- dim(.df)[2] # col of df
  pig_r_row <- dim(.pig_r)[1] # row of pigment ratio
  pig_r_col <- dim(.pig_r)[2] # col of df
  if (df_col != pig_r_col) {
    stop(paste(
      'normprod: matrix size error\nNot the same number of columns\ndf:',
      df_col,
      '\npig ratio:',
      pig_r_col)
    )
  }

  # ---- select column to normalize to ----
  # default = last column
  if (is.null(.col)) {
    # default to last col
    df_sel    <- ncol(.df)
    pig_r_sel <- ncol(.pig_r)
  } else {
    message(paste("Column selected for normalizing:", .col))
    df_sel    <- .col
    pig_r_sel <- .col
  }
  
  # ---- normalize df to inverse of last col ----
  df_norm_col    <- as.matrix(1 / (.df[,df_sel] + 1e-100))
  df_norm        <- .df * pracma::repmat(df_norm_col, 1, df_col)

  # ---- normalize pigment ratio to inverse of last col ----
  pig_r_norm_col <- as.matrix(.pig_r[,pig_r_sel])
  pig_r_norm     <- .pig_r * pracma::repmat(1 / (pig_r_norm_col + 1e-100), 1, df_col)
  
  # ---- normalize a * b to last col in df and pig ratio ----
  taxa_norm      <- .taxa_amt * pracma::repmat(df_norm_col, 1, pig_r_row) 
  cc             <- taxa_norm * pracma::repmat(t(pig_r_norm_col), df_row, 1)
  
  # ---- calc RMS ----
  res            <- df_norm - cc %*% pig_r_norm  # residuals
  rms            <- sqrt(mean(res^2))
  
  # ---- return final values ----
  results <- list(
                  ss = df_norm, 
                  cc = cc,
                  ff = pig_r_norm,
                  rms = rms
                  )
  
  return(results)
}
