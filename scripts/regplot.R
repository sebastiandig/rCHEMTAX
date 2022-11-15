regplot <- function(.df, .df_sd, .pig_r, .pig_r_sd, .info = NULL, verbose = TRUE,
                    fact_sd = NULL) {
################################################################################
#                                                                              # 
#          Plot the effect of weighting of f in positive matrix factors        #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# Runs `nnmatfactsd` with factors applied to the input `.pig_r_sd`. This records
# each iteration of the factors used and will have a resulting plot with the 
# factors used on the x axis and the weighted RMS on the y-axis. 
#
# ---- INPUTS: -----------
# df       = Matrix to be factored as a*b (best if a larger than b)
# df_sd    = Matrix of standard deviations for df
# pig_r    = Stabilizing value for pig_r, non zero locations, & initial value
# pig_r_sd = Matrix of standard deviations for pig_r
#
# ---- OUTPUTS: ----------
# results  = Results of the effects from the different factors on the starting 
#            pigment ratios matrix
#            $factors Vector of factors used to multiply pigment ratio matrix
#            $df      From `nnmatfactsd`, the factor multiplier and weighted RMS 
#                     for each factor
#            $plt     The plot `Weighted root mean square error vs standard 
#                     deviation x factor`
#            $logs    The data from each factor as output by `nnmatfactsd`
#
# ---- NOTES: ------------
# Original: 2010-02-21  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sat Jun 18 20:54:16 2022)
  
  # ---- load library ----
  library("tictoc")
  library("ggplot2")
  
  # set directory for saving 
  root <- rprojroot::find_rstudio_root_file()

  # source functions
  source(paste0(root,"/scripts/nnmatfactsd.R"))
  source(paste0(root,"/scripts/initstruct.R"))
  source(paste0(root,"/scripts/fancy_scientific.R"))
  
  # ---- create factors to test on pigment ratio matrix ----
  if (is.null(fact_sd)) {
    val     <- c(10, sqrt(10))
    fact_sd <-
      c(val, 0.1 * val, 0.01 * val, 0.02, 0.001 * val, 0.0001 * val)
  }
  
  n       <- length(fact_sd)
  rmsx    <- matrix(0, n, 1)
  rmsxwt  <- matrix(0, n, 1)
  
  # ---- initialize options ----
  if (is.null(.info)) .info  <- list()
  
  # ---- initialize list for logs ----
  logs <- list() 
  
  # ---- factor analysis with varying pig ratios multiplied by a factor  ----
  for (i in seq(fact_sd)) {
    
    cat(sprintf('\nFactor = %8.5f (%02d of %02d)\n', 
                    fact_sd[i], i, length(fact_sd)))
 
    tictoc::tic()

    # run factor analysis
    temp      <- nnmatfactsd(.df, 
                             .df_sd, 
                             .pig_r, 
                             .pig_r_sd * fact_sd[i], 
                             .info   = .info,
                             verbose = verbose)
 
    info      <- temp$info
    rmsx[i]   <- info$rmsx
    rmsxwt[i] <- info$rmsxwt
    
    # log each run
    fact_num         <- paste0("factor_", i)
    logs[[fact_num]] <- temp
    
    cat(sprintf('\n'))
    tictoc::toc()
    cat(sprintf('\n-----------------------------\n'))
    
  }    

  # ---- initialize plot info ----
  yticks_minor <- outer(1:10, 10^(-5:-1))
  xticks_minor <- outer(1:10, 10^(-4:1))
  
  df      <- data.frame(fact_sd = fact_sd, rmsxwt = rmsxwt)
  
  df_xmin <-  floor(log10(min(df$fact_sd)))
  
  # ---- print plot ----
  plt <- ggplot(df, aes(x = fact_sd, y = rmsxwt)) +
    
    geom_path() +
    geom_point(color = "red") +
    
    labs(
         title = 'Weighted root mean square error vs standard deviation x factor',
         x     = 'Factor on Pigment Ratio Standard Deviation',
         y     = 'Root mean square'
         ) +
    
    scale_x_log10(limits = c(10^(df_xmin), 100),
                  labels = fancy_scientific,
                  breaks = 10^(df_xmin:2),
                  minor_breaks = xticks_minor,
                  expand = expansion(add = c(0, 0.15))
                  ) +
    
    theme_bw() +
    theme(panel.grid.major.x = element_line(color = "gray"),
          panel.grid.minor.x = element_line(color = "gray")
          )
  
  if (verbose) {
    print(plt)
  }
  
  # ---- return final values ----
  results         <- list()
  results$factors <- fact_sd
  results$df      <- df
  results$plt     <- plt
  results$logs    <- logs
  
  return(results)
}