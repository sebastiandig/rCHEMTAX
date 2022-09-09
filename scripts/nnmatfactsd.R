nnmatfactsd <- function(.df,.df_sd,.pig_r_init,.pig_r_sd,.info=NULL, verbose = TRUE){
################################################################################
#                                                                              # 
#            Non negative matrix factors x=a*b with prior b0 & sd              #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# Can normalise matrices after function exits  eg  taxa_amt*inv(d), d*pig_r
# extends algorithm of Lee & Seung
#
# ---- INPUTS: -----------
# df         = Matrix to be factored as taxa_amt*pig_r (best if taxa_amt larger than pig_r)
# df_sd      = Matrix of standard deviations for x
# pig_r_init = Stabilizing value for pig_r, non zero locations, & initial value
#               (or initial pigment ratio)
# pig_r_sd   = Matrix of standard deviations for pig_r
# info       = Information/options for calculation (to override default values)
#              $maxitr   Maximum number of iterations
#              $convitr  Iterations between convergence tests
#              $printitr Count of iterations between printing
#              $conv     Converge test value for rms change
#              $inita    Initial taxa amount
#              $initb    Initial pigment ratio value
#              $maxaitr  Maximum for initial iterations for taxa_amt
#
# ---- OUTPUTS: ----------
# taxa_amt   = Left factor of x (i.e. amount of taxa per sample)
# pig_r      = Right factor of x
# info       = Initial info and
#              $itr         Iterations used
#              $conv        Converge level for error
#              $rms         Total rms error weighted for x & pig_r
#              $rms_df      Rms error in x unweighted
#              $rmdxwt      Weighted rms error in x
#              $rms_pig_r   Rms error in pig_r unweighted
#              $rmsbwt      Weighted rms error in pig_r
#         optional:
#           TODO: $plots
# ---- NOTES: ------------
# Original: 2010-03-12  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sun Jun 19 22:22:08 2022) 

  # ---- remove when done testing! ----
  # .df = s
  # .df_sd = ssd
  # .pig_r_init = f0
  # .pig_r_sd = fsd
  # .info = NULL
  
  
  library("pracma")
  
  root <- rprojroot::find_rstudio_root_file()
  
  source(paste0(root,"/scripts/initstruct.R"))
  source(paste0(root,"/scripts/amatfactsd.R"))
  
  # ---- options ----
  # create open list if no defaults are given for `info`  
  if (is.null(.info)) {
    .info <- list()
  }
  
  # set defaults options
  deflt <-
    list(
      maxitr   = 20000,
      convitr  = 500,
      printitr = 1000,
      conv     = 1e-6,
      initb    = .pig_r_init,
      maxaitr  = -1
    )
  
  # add/replace default options if not set
  info  <- initstruct(.info, deflt)
  
  # ensure iteration counts nest & convitr<=printitr<=maxitr
  convitr       <- info$convitr
  printitr      <- info$printitr
  maxitr        <- info$maxitr
  
  if (printitr > maxitr) {
    maxitr      <- ceiling(maxitr/convitr)*convitr
    # maxitr      <- pracma::ceil(maxitr/convitr)*convitr
    printitr    <- maxitr + 1 # if don't want to print values, could be different
  } else {
    printitr    <- ceiling(printitr/convitr)*convitr
    maxitr      <- ceiling(maxitr/printitr)*printitr
    # printitr    <- pracma::ceil(printitr/convitr)*convitr
    # maxitr      <- pracma::ceil(maxitr/printitr)*printitr
  }
  
  info$printitr <- printitr
  info$maxitr   <- maxitr
  
  # ---- extract row and column info ----
  # checks for same number of columns
  df_row           <- dim(.df)[1]  # row #
  df_col           <- dim(.df)[2]  # col #
  pig_r_init_row   <- dim(.pig_r_init)[1] # row #, # of classes to estimate 
  pig_r_init_col   <- dim(.pig_r_init)[2] # col #, # of pigments to adjust
  if (df_col != pig_r_init_col) stop('matfactsd: matrix size error\nnot the same size')
  
  # 
  pig_conc_tot     <- df_row * df_col            # total number of pigment conc for all samples
  tax_conc_tot     <- df_row * pig_r_init_row    # total number of taxa concentrations
  pig_adj_tot      <- sum(.pig_r_init != 0)      # number of pigments being adjusted
  n_err_terms      <- pig_conc_tot + pig_adj_tot # number of error terms (sum of adjusted pigments and pigment concentrations)
  
  # ---- add inverse-variance weighting variables ----
  df_sd            <- .df_sd + 1e-100     # create matrix of non-zeros
  pig_r_sd         <- .pig_r_sd  + 1e-100 # create matrix of non-zeros
  w2x              <- 1 / (df_sd^2)       # create weight by inverse squared std dev sample matrix
  w2b              <- 1 / (pig_r_sd^2)    # create weight by inverse squared std dev ratio matrix
  xw               <- .df * w2x           # multiply inverse-variance weighting x by samples
  b0w              <- .pig_r_init * w2b   # multiply inverse-variance weighting b by initial pigment ratio


  # ---- initial estimate before running factorization through iteration ----
  # initial estimate for b
  pig_r     <- info$initb
  
  # initial estimate for taxa_amt
  maxaitr   <- info$maxaitr
  
  # taxa_amt is:
  # 1: all contributions of taxa are random between 0 - 1
  # 2: all columns have the same contribution but different for each sample, based on 
  #    row sums of pigments
  if (exists("inita", info)) {
    taxa_amt       <- info$inita
  } else {
    taxa_amt       <-
      pracma::repmat(sqrt(as.matrix(apply(.df ^ 2, 1, sum)) / 
                            pig_r_init_row), 1, pig_r_init_row)
    if (maxaitr < 0) {
      maxaitr      <- 10
      info$maxaitr <- maxaitr
    }
  }

  # iterates starting taxa_amt, to get different starting contributions
  if (maxaitr > 0) {
    
    t1         <- xw %*%  t(pig_r)
    
    for (itr in 1:maxaitr) {
      taxa_amt <- 
        taxa_amt * t1 / (((taxa_amt %*% pig_r) * w2x) %*% t(pig_r) + 1e-100)
    }
    
    info$inita <- taxa_amt
  }

  # ---- initial rms of a and b estimates ----
  # initial pigments minus predicted pigments based on b (the pig ratio) and a (amount of pigment)
  df_res   <- .df - taxa_amt %*% pig_r 
  pig_r_res <- .pig_r_init - pig_r
  
  # error for a and b, from initial divided by std dev
  res.std   <- rbind(df_res / df_sd, pig_r_res / pig_r_sd ) # this is standardized residuals
  rms_prev  <- sqrt(sum(res.std^2) / n_err_terms) # initial root mean square 
  conv      <- info$conv
  
  # ---- RMS of initial pigments to predicted pigments ----
  rms_df   <- sqrt(sum(df_res^2)/pig_conc_tot)
  rmsxwt    <- sqrt(sum((df_res/df_sd)^2/pig_conc_tot))
  
  # alternate:
  # rms_df <- Metrics::rmse(df, taxa_amt %*% pig_r)
  # rmsxwt  <- Metrics::rmse(df/df_sd, (taxa_amt %*% pig_r)/df_sd)
  
  # residuals of initial pigment ratio to current pigment ratio
  rms_pig_r <- sqrt(sum(pig_r_res^2) / pig_adj_tot) # rms of pigment ratio
  rmsbwt    <- sqrt(sum((pig_r_res / pig_r_sd )^2) / pig_adj_tot)
    
  
  # ---- print heading & initial values ----
  if (printitr <= maxitr & verbose) { 
    cat(sprintf('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:     dRMS:   dTaxa RMS:   dPig RMS:\n'))
    cat(sprintf("%5i%#8.3g%#10.3g%#21.4g%#19.3g\n", 0,rms_df,rmsxwt,rms_pig_r,rmsbwt))
  }
  
  # ---- initialize factorization ----
  # set previous iteration values to check against
  taxa_amt_prev <- taxa_amt
  pig_r_prev    <- pig_r

  # initialize log
  # TODO: be an option, should add columns for each species
  logs <-
    data.frame(
      itr          = 0,
      rms_df      = rms_df,
      rmsxwt       = rmsxwt,
      rms_pig_r    = rms_pig_r,
      rmsbwt       = rmsbwt,
      rms_chg      = NA,
      taxa_amt_chg = NA,
      pig_r_chg    = NA
    )

  # ---- main factorization loop ----
  for (itr in seq(maxitr)) {
    # update a & b
    taxa_amt <- taxa_amt * (xw %*% t(pig_r)) / (((taxa_amt %*% pig_r) * w2x) %*% t(pig_r) + 1e-100)
    pig_r    <-
      pig_r * ((t(taxa_amt) %*% xw) + b0w) / (t(taxa_amt) %*% ((taxa_amt %*% pig_r) * w2x) + pig_r * w2b + 1e-100)
    
    # check convergence occasionally
    if (itr %% convitr == 0) {
      
      # change in residuals from previous # of iterations 
      taxa_amt_chg  <- sqrt(sum((taxa_amt_prev - taxa_amt)^2) / tax_conc_tot) / convitr
      pig_r_chg     <- sqrt(sum((pig_r_prev - pig_r)^2) / pig_adj_tot) / convitr
      
      # calc residuals divided by std. dev.
      df_res       <- .df - taxa_amt %*% pig_r 
      pig_r_res     <- .pig_r_init - pig_r
      res.std       <- rbind(df_res / df_sd, pig_r_res / pig_r_sd )
      
      # calc rms and rms change from previous iterations divided by # of iterations since last time (default = 500)
      rms           <- sqrt(sum(res.std^2) / n_err_terms)
      rms_chg       <- (rms_prev - rms ) / convitr # change in rms from previous to current, where divide by nuber of convitr
      
      # non-negative least linear squares to calc matrix `a` in min||x - a %*% b||
      taxa_amt      <- amatfactsd(.df,df_sd,pig_r)
      taxa_amt_prev <- taxa_amt
      pig_r_prev    <- pig_r
      rms_prev      <- rms
      
      # when the change in RMS is < conv, then breaks
      conv_end      <- abs(rms_chg) < conv

      # check print occasionally
      if ((itr %% printitr) == 0 || conv_end || (itr %% maxitr) == 0  ) {
        df_res   <- .df - taxa_amt %*% pig_r
        rms_df   <- sqrt(sum(df_res^2) / pig_conc_tot)
        rmsxwt    <- sqrt(sum((df_res / df_sd)^2 / pig_conc_tot))
        pig_r_res <- .pig_r_init - pig_r
        rms_pig_r <- sqrt(sum(pig_r_res^2) / pig_adj_tot)
        rmsbwt    <- sqrt(sum((pig_r_res / pig_r_sd )^2) / pig_adj_tot)
        
        logs <- rbind(logs, cbind(itr,rms_df,rmsxwt,rms_pig_r,rmsbwt,rms_chg,taxa_amt_chg,pig_r_chg))
        
        if (printitr <= maxitr & verbose) {
          cat(sprintf('%5i%#8.3g%#10.3g%#21.4g%#19.2f%#10.2e%13.2e%#12.2e\n', 
                       itr,rms_df,rmsxwt,rms_pig_r,rmsbwt,rms_chg,taxa_amt_chg,pig_r_chg) )
          cat(sprintf("%5i%#8.3g%#10.3g%#21.4g%#19.3g\n", 0,rms_df,rmsxwt,rms_pig_r,rmsbwt))
        }
      }
      # pracma::fprintf( "%#19.5f", 12. )
      # check for convergence
      if (conv_end) {
        break
      }
    }
  }
  
  # ---- return final values ----
  info$itr       <- itr       # final iteration
  info$conv      <- rms_chg   # convergence value
  info$rms       <- rms       # total rms
  info$rms_df   <- rms_df     # df rms
  info$rmsxwt    <- rmsxwt 
  info$rms_pig_r <- rms_pig_r # pigment ratio rmc
  info$rmsbwt    <- rmsbwt
  
  info$logs <- logs
  # TODO: make option for plots
  # if (isTRUE(plots)) {
    # info$plots <-  # maybe inlcude this plot at the end?
    #   ggplot(logs, aes(itr, rms_chg)) + 
    #   geom_line(color = "black") + 
    #   geom_point() + 
    #   geom_hline(yintercept = tes$info$conv, color = "red") + 
    #   scale_y_log10()
  # }
  
  
  result <- list(a = taxa_amt, b = pig_r, info = info)
  
  return(result)
  
  # ---- end ----  
}