randstart <- function(.df,.df_sd,.pig_r,.pig_r_sd, .info = NULL,.nrep = 10, verbose = TRUE) {
################################################################################
#                                                                              # 
#   Plot the effect of random starts for s% effect of random start locations   #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# 
# ---- INPUTS: -----------
# df       = Matrix to be factored as a*b (best if a larger than b)
# df_sd    = Matrix of standard deviations for x
# pig_r    = Stabilizing value for b, non zero locations, & initial value
# pig_r_sd = Matrix of standard deviations for b
#
# ---- OUTPUTS: ----------
# TODO: edit this
#
# ---- NOTES: ------------
# Original: 2010-02-26  Matlab7  W.Whiten
#  
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sun Jun 19 14:01:35 2022)  

  # ---- remove after testing! ----
  # .nrep <- 10
  
  # ---- load library ----
  library("tictoc")
  library("pracma")
  library("ggplot2")
  
  # set directory for saving 
  root <- rprojroot::find_rstudio_root_file()
  
  # references another function
  source(paste0(root,"/scripts/nnmatfactsd.R"))
  source(paste0(root,"/scripts/initstruct.R"))
 
  # ---- initialize matrices for random starting points ----
  df_row    <- dim(.df)[1]             # row # .df
  pig_r_col <- dim(.pig_r)[1]          # row # .pig_r 
  indx      <- which(.pig_r > 0)       # index of non-zero pigments
  n_pig     <- length(indx)            # number of pigments
  
  # initialize pigment ratio and taxa contributions matrix for each random start
  pig_rep   <- matrix(0, .nrep, n_pig) 
  taxa_amt_rep <- matrix(0, .nrep, df_row*pig_r_col)
  
  # initialize list for logs
  logs <- list() 
  
  if (is.null(.info)) {
    info_init <- list()
  } else {
    info_init <- .info
  }
  
  # ---- random starting points loop ----
  # starts clock to see time passed
  tictoc::tic.clearlog()
  start <- tictoc::tic()
  
  for (i in seq(.nrep)) {
    tictoc::tic()
    
    # ---- initialize options ----
    # TODO: make sure original settings do not get overwritten with each iteration
    deflt <- list(
                    inita    = pracma::rand(df_row, pig_r_col),
                    maxitr   = 30000,
                    printitr = 1e12,
                    conv     = 1e-10
                  )
    
  
    # add/replace default options if not set
    .info  <- initstruct(info_init, deflt)
    
    # ---- run factor analysis ----
    cat(sprintf('\nRandom Start Number: %02d of %02d\n', i, .nrep))
    # pracma::fprintf('\nRandom Start Number: %02d of %02d\n', i, .nrep)

    temp          <- nnmatfactsd(.df,
                                 .df_sd,
                                 .pig_r,
                                 .pig_r_sd,
                                 .info   = .info,
                                 verbose = verbose
                                 )
    taxa_amt_temp <- temp$a
    pig_r_temp    <- temp$b
    info          <- temp$info

    # log each run
    rep_num          <- paste0("rep_num_", i)
    logs[[rep_num]]  <- temp

    # log taxa contribution and pigment ratio out
    pig_rep[i,]      <- t(pig_r_temp[indx])
    taxa_amt_rep[i,] <- t(taxa_amt_temp)

    end <- (tictoc::toc(quiet = T))$toc - start

    # when printitr is above maxitr, will only print last iterations
    if (!(info$printitr <= info$maxitr) & verbose) {
      cat(sprintf('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:     dRMS:   dTaxa RMS:   dPig RMS:\n'))
      cat(sprintf('%5i%#8.3g%#10.3g%#21.3g%#19.3g%#10.2e%13.2e%#12.2e\n',
                  info$itr,info$rms_pig,info$rmsxwt,info$rms_pig_r,info$rmsbwt,
                  info$logs$rms_chg[2],info$logs$taxa_amt_chg[2],info$logs$pig_r_chg[2]))
      
      # pracma::fprintf('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:     dRMS:   dTaxa RMS:   dPig RMS:\n')
      # pracma::fprintf('%5i%#8.3g%#10.3g%#21.3g%#19.3g%#10.2e%13.2e%#12.2e\n',
      #                 info$itr,info$rms_pig,info$rmsxwt,info$rms_pig_r,info$rmsbwt,
      #                 info$logs$rms_chg[2],info$logs$taxa_amt_chg[2],info$logs$pig_r_chg[2])
    }
    
    cat(sprintf('\n-----------------------------\n'))
    cat(sprintf('\nTime Elapsed: %.2f (s)\nConverged at: %.4e\n',
                end, info$conv))
    cat(sprintf('\n-----------------------------\n'))
    # pracma::fprintf('\n-----------------------------\n')
    # pracma::fprintf('\nTime Elapsed: %.2f (s)\nConverged at: %.4e\n',
    #                 end, info$conv)
    # pracma::fprintf('\n-----------------------------\n')
    
  }
  
  # ---- display range of values for each pigment and taxa contributions ----
  pig_rep_range    <-
    pig_rep - pracma::repmat(apply(pig_rep, 2, mean, na.rm = T), .nrep, 1)
  
  df_pig_rep_range <-
    taxa_amt_rep - pracma::repmat(apply(taxa_amt_rep, 2, mean, na.rm = T), .nrep, 1)
  
  if (verbose) {
    cat(sprintf('\nRange of variation in Pigment Ratios and Taxa Contribution:\n'))
    cat(sprintf('\n-----------------------------\n'))
    cat(sprintf('\nPigment Ratio:\n%14s%7.3g\n%14s%7.3g\n',
                "Max:",min(pig_rep_range, na.rm = T), "Min:", max(pig_rep_range, na.rm = T)))
    cat(sprintf('\nTaxa Contributon:\n%14s%7.3g\n%14s%7.3g\n',
                "Max:",min(df_pig_rep_range, na.rm = T), "Min:",max(df_pig_rep_range, na.rm = T)))
    cat(sprintf('\n-----------------------------\n'))
    
    # pracma::fprintf('\nRange of variation in Pigment Ratios and Taxa Contribution:\n')
    # pracma::fprintf('\n-----------------------------\n')
    # pracma::fprintf('\nPigment Ratio:\n%14s%7.3g\n%14s%7.3g\n',
    #                 "Max:",min(pig_rep_range, na.rm = T), "Min:", max(pig_rep_range, na.rm = T))
    # pracma::fprintf('\nTaxa Contributon:\n%14s%7.3g\n%14s%7.3g\n',
    #                 "Max:",min(df_pig_rep_range, na.rm = T), "Min:",max(df_pig_rep_range, na.rm = T))
    # pracma::fprintf('\n-----------------------------\n')
  }
  
  # ---- initialize plotting variables ----
  source(paste0(root,"/scripts/fancy_scientific.R"))
  yticks_minor = outer((1:10),(10^(-5:-1)))
  xticks_minor = outer((1:10),(10^(0:1)))
  
  
  # figure
  # loglog(mean(tf),std(tf)./mean(tf),'o')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
  # pigment ratios
  df <- data.frame(x = apply(pig_rep, 2, mean, na.rm = TRUE),
                   y = apply(pig_rep, 2, function(x) 
                     sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
                   )
  
  df_xmin <- floor(log10(min(df$x)))
  df_ymin <- floor(log10(min(df$y)))
  
  # taxa contribution
  # might have issue with transpose?
  df2 <- data.frame(x = apply(taxa_amt_rep + 1e-100, 2, mean, na.rm = TRUE),
                    y = apply(taxa_amt_rep + 1e-100, 2, function(x)
                      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  )
  
  df2_xmin <- floor(log10(min(df2$x)))
  df2_ymin <- floor(log10(min(df2$y)))
  
  # ---- create plots ----
  # pigment ratios
  plt1 <- ggplot(data = df, aes(x = x, y = y)) + 
  geom_point(colour = "red") +
    # scale_x_log10(limits = c(10^df_xmin, 100),
    #               labels = fancy_scientific,
    #               breaks = 10^(df_xmin:2),
    #               minor_breaks = xticks_minor,
    #               expand = c(0,0)) +

    # scale_y_log10(limits = c(10^(df_ymin), 100),
    #               labels = fancy_scientific,
    #               breaks = c(10^(df_ymin:2)),
    #               minor_breaks = yticks_minor,
    #               expand = c(0,0)) +
    labs(
         title = "Variation in Pigment Ratio Coefficients",
         x     = "Mean Coefficient Value",
         y     = "Coefficient of Variation (Standard Deviation / Mean)") +
    theme_bw()
  
  
  # figure
  # loglog(mean(tc)',(std(tc)./mean(tc))','.')
  # title('Variation in c coefficients')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
  # taxa contribution
  plt2 <- ggplot() +
    # original had bubble dots - shape?
    geom_point(data = df2,
               aes(x = x,
                   y = y),
               colour = "red") +
      # scale_x_log10(limits = c(10^df2_xmin, 100),
      #               labels = fancy_scientific,
      #               breaks = 10^(df2_xmin:2),
      #               minor_breaks = xticks_minor,
      #               expand = c(0,0)) +
      #
      # scale_y_log10(limits = c(10^(df2_ymin), 100),
      #               labels = fancy_scientific,
      #               breaks = c(10^(df2_ymin:2)),
      #               minor_breaks = yticks_minor,
      #               expand = c(0,0)) +
    labs(title = "Variation in Taxa Contribution",
         x     = "Mean Coefficient Value",
         y     = "Standard Deviation / Mean") +
    theme_bw()
  
  if (verbose) {
    print(plt1)
    print(plt2)
  }
  
  
  # ---- return final values ----
  results            <- list()
  results$replicates <- .nrep
  results$rep_df     <- list(pig_rep = pig_rep, taxa_amt_rep = taxa_amt_rep)
  results$dfs        <- list(pig_r_avg = df,taxa_r_avg = df2)
  results$plt        <- list(pig_plt = plt1, taxa_plt = plt2)
  results$logs       <- logs
  
  return(results)
  }
