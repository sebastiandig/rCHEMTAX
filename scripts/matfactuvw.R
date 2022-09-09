matfactuvw <- function(.df, u = NULL, v = NULL, .pig_r, w = NULL, wb = 1, .info = NULL, verbose = TRUE) {
################################################################################
#                                                                              # 
#               Non neg matrix factor with factored sd values                  #
#              Fast method in Matlab with approx standard deviations           #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# 
# ---- INPUTS: -----------
# df       = Matrix to be factored as a*.pig_r
# u       = Row sd factor for df (sd is u'*v)
# v       = Column sd factor for df
# pig_r = Stabilising value for b, non zero locations, & initial value
# w       = Row sd factor for for (sd is w'*v)
# wb      = Scalar weight factor for b errors (short call sequence only)
# info    = Information for calculation (to override default values)
#           $maxitr  Maximum number of iterations
#           $printitr Count of iterations between printing
#           $conva  Converge value for a
#           $convb  Converge value for b
#           $conve  Converge value for rms change
#           $inita  Initial a value
#           $initb  Initial b value
#
# ---- OUTPUTS: ----------
# a       = Left factor of df
# b       = Right factor of df
# info    = Initial info and
#           $itr  Iterations used
#           $conva Converge level for a
#           $convb Converge level for b
#           $conve Converge level for error
#
# ---- NOTES: ------------
# Original: 2010-02-21  Matlab7  W.Whiten
# 
# ---- REFERENCES(s): -----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Mon Jun 20 21:26:19 2022)
  
  # ---- TODO: remove after testing ----
  # x = s_norm
  # u = NULL
  # v = NULL
  # b = f_norm
  # w = NULL
  # wb = 1
  # .info = NULL
  # verbose = T
  # 
  # .df = s_norm
  # .pig_r = f_norm
  # 
  # ---- set directory ----
  root <- rprojroot::find_rstudio_root_file()
  source(paste0(root,"/scripts/initstruct.R"))
  
  # ---- extract row and column info ----
  # checks for same number of columns
  df_row           <- dim(.df)[1]  # row #
  df_col           <- dim(.df)[2]  # col #
  pig_r_init_row   <- dim(.pig_r)[1] # row #, # of classes to estimate 
  pig_r_init_col   <- dim(.pig_r)[2] # col #, # of pigments to adjust
  if (df_col != pig_r_init_col) stop('matfactsd: matrix size error\nnot the same size')
  
  # number of pigments
  pig_tot <- sum(.pig_r != 0)
  
  # ---- set up variable if not given ----
  if (is.null(u)) {
    # u=sqrt(sum(x.^2,2)/size(x,2));
    u <- sqrt(apply(.df^2, 1, sum, na.rm = T) / df_col)
    u <-  matrix(u)
    # u <-  (apply(.df, 1, sd, na.rm = T)) # <- this seems more accurate?
  }
  
  if (is.null(v)) {
    # v=0.01*sqrt(sum(x.^2,1)/size(x,1))+0.01;
    v <- 0.01 * sqrt(apply(.df^2, 2, sum, na.rm = T) / df_row) + 0.01
    v <- matrix(v)
  }
  
  if (is.null(w)) {
    w      <- wb * sqrt(apply(.pig_r^2, 1, sum, na.rm = T) / pig_r_init_col)
    w_diag <- diag(1 / w)
    w      <-  matrix(w)
  }

  # ---- options ----
  # create open list if no defaults are given for `info`  
  if (is.null(.info)) {
    .info <- list()
  }
  
  # set defaults options
  deflt <- list(
    maxitr   = 1000,
    # convitr  = 500, # TODO: include this
    printitr = 100,
    conve    = 1e-6,
    conva    = 1e-6,
    convb    = 1e-6,
    inita    = matrix(1, df_row, pig_r_init_row) / pig_r_init_row,
    initb    = .pig_r
  )
  
  # add/replace default options if not set
  info  <- initstruct(.info, deflt)
  
  maxitr <- info$maxitr
  printitr <- info$printitr
  # convitr <- info$conitr
  
  
  # ---- initial rms and estimates of a and b ----
  # initial a
  taxa_amt <- info$inita / pracma::repmat(u,1,pig_r_init_row)
  
  # initial b
  pig_r <- info$initb / pracma::repmat(t(v),pig_r_init_row,1)
  pig_r <- pig_r / pracma::repmat(as.matrix(v[length(v)] * pig_r[,ncol(pig_r)]),1,df_col)
  
  # data divided by std deviation of col and row
  df_uv    <- .df / (u %*% t(v))
  df_pig_r <- rbind(df_uv, f_norm / w %*% t(v))
  
  # calc initial residuals and rms
  res <- df_pig_r - rbind(taxa_amt, w_diag) %*% pig_r
  rms <- sqrt( sum(res^2) / length(res))
  
  
  # ---- RMS of initial pigments to predicted pigments ----
  # calc residuals
  df_uv_res <- df_uv - taxa_amt %*% pig_r
  pig_r_res <- .pig_r / (w %*%  t(v)) - w_diag %*% pig_r
  
  taxa_t  <- taxa_amt * pracma::repmat(u, 1, pig_r_init_row)
  pig_r_t <- pig_r * pracma::repmat(t(v), pig_r_init_row, 1)
  
  df_res <- .df - taxa_t %*% pig_r_t 
  
  # calc rms for df and pigment ratio
  rms_df_uv <- sqrt(mean(df_uv_res^2))
  rms_pig_r <- sqrt(sum(pig_r_res^2) / pig_tot)
  rms_df   <-  sqrt(mean(df_res^2))

  # ---- print heading & initial values ----
  if (info$printitr < 1e6 & verbose) {
    # cat(sprintf('   itr    rmsx      rmsxwt     rmsb     drms       daa       dbb\n'))
    # cat(sprintf('%6i%#11.3g%#11.3g%#11.3g\n',0,rms_df,rms_df_uv,rms_pig_r))
    
    cat(sprintf('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:      dRMS:   dTaxa RMS:   dPig RMS:\n'))
    cat(sprintf("%5i%#8.3g%#10.4g%#21.4g\n", 0,rms_df,rms_df_uv,rms_pig_r))
    }
  
  # ---- initialize log ----
  # TODO: be an option, should add columns for each species
  logs <-
    data.frame(
      itr          = 0,
      rms_df       = rms_df,
      rms_df_uv    = rms_df_uv,
      rms_pig_r    = rms_pig_r,
      # rmsbwt     = rmsbwt,
      rms_chg      = NA,
      taxa_amt_chg = NA,
      pig_r_chg    = NA
    )
  
  # ---- main factorization loop ----
  # aa1=aa.*(xuv*bb')./(aa*(bb*bb'));
  # %aa1=aa.*(1+((xuv-aa*bb)*bb')./(aa*(bb*bb')));
  
  for (itr in seq(maxitr)) {
    taxa_amt_cur <- taxa_amt * (df_uv %*% t(pig_r)) / (taxa_amt %*% (pig_r %*% t(pig_r)))
  
    # bb1=bb.*(aw'*xb)./(aw'*(aw*bb));
    # %bb1=bb.*(1+(aw'*(xb-aw*bb))./(aw'*(aw*bb)));
    # %bb1=bb1./repmat(v(end).*bb(:,end),1,np);
    
    taxa_amt_wt  <- rbind(taxa_amt_cur, w_diag)
    pig_r_cur    <-
      pig_r * (t(taxa_amt_wt) %*% df_pig_r) / (t(taxa_amt_wt) %*% (taxa_amt_wt %*% pig_r))
    
    # ---- calc change in rms ----
    taxa_amt_chg <- sqrt(mean((taxa_amt_cur - taxa_amt)^2))
    pig_r_chg    <- sqrt(sum((pig_r_cur - pig_r)^2) / pig_tot)
    res          <- df_pig_r - taxa_amt_wt %*% pig_r_cur
    rms_prev     <- sqrt(mean(res^2))
    rms_chg      <- rms - rms_prev 
    
    # ---- set prev rms ----
    taxa_amt  <- taxa_amt_cur
    pig_r  <- pig_r_cur
    rms <- rms_prev
    
    # check convergence occasionally
    # if (itr %% convitr == 0) { # TODO: include convitr
    if (itr %% printitr == 0) {
      
      # calc residuals
      df_uv_res <- df_uv - taxa_amt %*% pig_r
      pig_r_res <- .pig_r / (w %*%  t(v)) - w_diag %*% pig_r
      
      taxa_t  <- taxa_amt * pracma::repmat(u, 1, pig_r_init_row)
      pig_r_t <- pig_r * pracma::repmat(t(v), pig_r_init_row, 1)
      df_res <- .df - taxa_t %*% pig_r_t 
      
      # calc rms for df and pigment ratios
      rms_df_uv <- sqrt(mean(df_uv_res^2))
      rms_pig_r <- sqrt(sum(pig_r_res^2) / pig_tot)
      rms_df   <-  sqrt(mean(df_res^2))
      
      logs <- rbind(logs, cbind(itr,rms_df,rms_df_uv,rms_pig_r,rms_chg,taxa_amt_chg,pig_r_chg))
      
      if (printitr <= maxitr & verbose) {
        cat(sprintf('%5i%#8.3g%#10.4g%#21.4g%#11.2e%13.2e%#12.2e\n', 
                    itr,rms_df,rms_df_uv,rms_pig_r,rms_chg,taxa_amt_chg,pig_r_chg))
      }
      
    }
    
    if (taxa_amt_chg < info$conva | pig_r_chg < info$convb | rms_chg < info$conve) {
      if (itr %% printitr == 0) {
        break
      } else {
        logs <- rbind(logs, cbind(itr,rms_df,rms_df_uv,rms_pig_r,rms_chg,taxa_amt_chg,pig_r_chg))
        break
      }
    }
  }
  
  # ---- return final values ----
  taxa_final <- taxa_amt * pracma::repmat(u, 1, pig_r_init_row)
  pig_r_fin <- pig_r * pracma::repmat(t(v), pig_r_init_row, 1)
  chlor_norm <- matrix(pig_r_fin[, ncol(pig_r_fin)])
  pig_r_fin <- pig_r_fin / pracma::repmat(chlor_norm, 1, df_col)
  taxa_norm <- taxa_final *  pracma::repmat(t(chlor_norm), df_row, 1)
  
  info$itr       <- itr
  info$conva     <- taxa_amt_chg
  info$convb     <- pig_r_chg
  info$conve     <- rms_chg
  info$rms_df    <- rms_df
  info$rms_df    <- rms_df
  info$rms_df_uv <- rms_df_uv
  info$rms_pig_r <- rms_pig_r

  
  results <- list(a = taxa_norm, b = .pig_r, chlor_norm = chlor_norm, info, logs = logs)     
  
  return(results)    
  # ---- end ----
}