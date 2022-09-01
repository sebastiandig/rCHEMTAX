regplot <- function(.df,.df_sd,.pig_r,.pig_r_sd) {
################################################################################
#                                                                              # 
#          Plot the effect of weighting of f in positive matrix factors        #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# Ploting function
#
# ---- INPUTS: -----------
# df       = Matrix to be factored as a*b (best if a larger than b)
# df_sd    = Matrix of standard deviations for df
# pig_r    = Stabilizing value for pig_r, non zero locations, & initial value
# pig_r_sd = Matrix of standard deviations for pig_r
#
# ---- OUTPUTS: ----------
# Plot
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
  
  source(paste0(root, "/scripts/nnmatfactsd.R"))
  
  # ---- create factors to test on pigment ratio matrix ----
  val    <- c(10, sqrt(10))
  fsdx   <- c(val, 0.1*val, 0.01*val, 0.02, 0.001*val, 0.0001*val)
  n      <- length(fsdx)
  rmsx   <- matrix(0, n, 1)
  rmsxwt <- matrix(0, n, 1)
  
  logs <- list()
  
  # ---- factor analysis with varying pig ratios multipled by a factor  ----
  for (i in seq(fsdx)) {
    
    
    pracma::fprintf('\nFactor = %8.5f (%02d of %02d)\n\n', fsdx[i], i, length(fsdx))
    
    tictoc::tic()
    temp      <- nnmatfactsd(.df,.df_sd,.pig_r,fsdx[i]*.pig_r_sd)
    
    info      <- temp$info
    rmsx[i]   <- info$rmsx
    rmsxwt[i] <- info$rmsxwt
    
    
    fact_num         <- paste0("factor_", i)
    logs[[fact_num]] <- temp
    
    pracma::fprintf('\n')
    tictoc::toc()
    pracma::fprintf('\n-----------------------------\n')
  }    

  # ---- initialize plot info ----
  source(paste0(root,"/scripts/fancy_scientific.R"))
  yticks = outer(1:10, 10^(-5:-1))
  xticks_minor = outer(1:10, 10^(0:1))
  
  df <- data.frame(fsdx = fsdx,rmsxwt = rmsxwt)
  
  df_xmin <-  floor(log10(min(df$fsdx)))
  
  # ---- print plot ----
  plt <- ggplot(df, aes(x = fsdx, y = rmsxwt)) +
    
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
  
  print(plt)
  
  results         <- list()
  results$factors <- fsdx
  results$df      <- df
  results$plt     <- plt
  results$logs    <- logs
  
  return(results)
}