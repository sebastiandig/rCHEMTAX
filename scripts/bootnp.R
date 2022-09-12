bootnp <- function(.df,.df_sd,.pig_r,.pig_r_sd, .info =NULL,.nrep = 10, 
                   verbose = TRUE, .pigm = NULL, .taxa = NULL) {
################################################################################
#                                                                              # 
#           Plot the effect of non parametric bootstrap of s                   #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
#
# ---- INPUTS: -----------
# .df   =  Matrix to be factored as a*b (best if a larger than b)
# .df_sd  = Matrix of standard deviations for x
# .pig_r  = Stabilising value for b, non zero locations, & initial value
# .pig_r_sd = Matrix of standard deviations for b
#
# TODO: need to update
# 
# ---- OUTPUTS: ----------
# TODO: need to update
#
# ---- NOTES: -------------
# Original: 2010-02-27  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Fri Jun 10 12:01:31 2022)
  
  # ---- TODO: remove after testing ----
  # TODO: graphing part and figure out cc1
  .df = s
  .df_sd = ssd
  .pig_r = f0
  .pig_r = fsd
  .info = NULL
  .info = list(printitr = 1000)
  .nrep = 2
  verbose = T
  .pigm = NULL
  .taxa = NULL
  
  # ---- load library ----
  library("tictoc")
  library("pracma")
  library("ggplot2")
  
  # set directory for saving 
  root <- rprojroot::find_rstudio_root_file()

  # source functions
  source(paste0(root,"/scripts/nnmatfactsd.R"))
  source(paste0(root,"/scripts/fancy_scientific.R"))
  source(paste0(root,"/scripts/initstruct.R"))

  # ---- initialize matrices ----
  df_row       <- dim(.df)[1]    # row # .df
  df_col       <- dim(.df)[2]    # col # of .df
  pig_r_row    <- dim(.pig_r)[1] # row of .pig_r
  pig_r_col    <- dim(.pig_r)[2] # col of .pig_r
  
  indx         <- which(.pig_r > 0)  # index of non-zero pigments
  n_pig        <- length(indx)       # number of non-zero pigments
  pig_r_rep    <- matrix(0, .nrep, n_pig)
  taxa_amt_rep <- matrix(0, .nrep, df_row*pig_r_row)

  # initialize empty matrix
  # cc1=zeros(df_row,pig_r_row); ????
  cc1  <- matrix(0,df_row,pig_r_row) 

  # ---- options ----
  # create open list if no defaults are given for `info`  
  if (is.null(.info)) .info  <- list()
  
  # set defaults options
  deflt    <- list(printitr = 1e6)
  
  # add/replace default options if not set
  info     <- initstruct(.info, deflt)
  
  # ---- initialize list for logs ----
  logs <- list() 
  
  # ---- start loop ----
  tictoc::tic.clearlog()
  start <- tictoc::tic()
  
  for (i in seq(.nrep)) {
    tictoc::tic()
    
    cat(sprintf('\nRandom Start Number: %02d of %02d\n', i, .nrep))
    
    # select subset of data
    ind <- ceiling(pracma::rand(df_row,1) * df_row)

    temp   <-
      nnmatfactsd(.df[ind, ], .df_sd[ind, ], .pig_r, .pig_r_sd, .info = info)
    taxa_amt_temp <- temp$a
    pig_r_temp    <- temp$b
    info          <- temp$info
    # only weirdness is cc1
    cc1[,] <- NA
    cc1[ind,] <- taxa_amt_temp
    
    # log taxa contribution and pigment ratio out
    pig_r_rep[i,]    <- t(pig_r_temp[indx])
    taxa_amt_rep[i,] <- t(cc1)
    # taxa_amt_rep[i,] <- t(taxa_amt_temp)
    
    # log each run
    rep_num          <- paste0("rep_num_", i)
    logs[[rep_num]]  <- temp
    
    end <- (tictoc::toc(quiet = T))$toc - start
    
    # when printitr is above maxitr, will only print last iterations
    if (!(info$printitr <= info$maxitr) & verbose) {
      cat('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:    ',
          'dRMS:   dTaxa RMS:   dPig RMS:',
          sprintf('\n%5i%#8.3g%#10.3g%#21.3g%#19.3g%#10.2e%13.2e%#12.2e\n',
                  info$itr,info$rms_pig,info$rmsxwt,info$rms_pig_r,info$rmsbwt,
                  info$logs$rms_chg[2],info$logs$taxa_amt_chg[2],
                  info$logs$pig_r_chg[2]))
    }
    
    cat(sprintf('\n-----------------------------
      \nTime Elapsed: %.2f (s)\nConverged at: %.4e
      \n-----------------------------\n',
      end, info$conv)
    )
  }

  # ---- display range of values for each pigment and taxa contributions ----
  pig_rep_range <-
    pig_r_rep - pracma::repmat(apply(pig_r_rep, 2, mean, na.rm = T), .nrep, 1)
  taxa_amt_rep_range <-
    taxa_amt_rep - pracma::repmat(apply(taxa_amt_rep, 2, mean, na.rm = T), .nrep, 1)
  
  if (verbose) {
    
    cat('\nRange of variation in Pigment Ratios and Taxa Contribution:
      \n-----------------------------\n',
      sprintf('\nPigment Ratio:\n%14s%7.3f\n%14s%7.3f\n',
              "Max:", min(pig_rep_range, na.rm = T), 
              "Min:", max(pig_rep_range, na.rm = T)),
      sprintf('\nTaxa Contributon:\n%14s%7.3f\n%14s%7.3f\n',
              "Max:", min(taxa_amt_rep_range, na.rm = T), 
              "Min:", max(taxa_amt_rep_range, na.rm = T)),
      '\n-----------------------------\n')
    
  }
  
  # ---- calc avg and sd for final matrices ----
  # pig_r_rep
  pig_r_fin <- 
    data.frame(
      x      = apply(pig_r_rep, 2, mean, na.rm = TRUE),
     std_dev = apply(pig_r_rep, 2, sd, na.rm = TRUE),
     y       = apply(pig_r_rep, 2, function(x) 
                     sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
  )
  
  
  # I don't think I need this chunk, but maybe
  # mtc=zeros(ns*nt,1);
  # stc=zeros(ns*nt,1);
  # for i=1:ns*nt
  #   t=tc(~isnan(tc(:,i)),i);
  #   mtc(i)=mean(t);
  #   stc(i)=std(t);
  # end
  
  taxa_amt_fin <- 
    data.frame(
      x       = apply(taxa_amt_rep, 2, mean, na.rm = TRUE),
      std_dev = apply(taxa_amt_rep, 2, sd, na.rm = TRUE),
      y       = apply(df, 2, function(x)
                      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  )
  
  # figure
  # loglog(mean(tf),std(tf)./mean(tf),'o')
  # title('Non parametric bootstrap: Variation in f coefficients')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
  # ---- initialize plot info ----
  yticks_minor = outer(1:10, 10^(-5:-1))
  xticks_minor = outer(1:10, 10^(0:1))

  
  
  plt_pig_r <- ggplot() +
    # original had bubble dots - shape?
    geom_point(data = pig_r_fin,
               aes(x = x,
                   y = y),
               colour = "red") +
    # scale_x_log10(limits = c(1, 100),
    #               labels = fancy_scientific,
    #               minor_breaks = xticks) +
    # scale_y_log10(limits = c(1e-5, 1),
    #               labels = fancy_scientific,
    #               minor_breaks = yticks) +
    labs(
      title = paste("Non Parametric Bootstrap:",
                    "Variation in Pigment Ratio coefficients"),
      x     = "Mean coefficient value",
      y     = "Coefficient of Variation"
    ) +
    theme_bw()

# figure
# %loglog(mtc',stc,'.')
# loglog(mtc,stc./mtc,'.')
# title('Non parametric bootstrap: Variation in c coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')

  # taxa_amt_rep
  plt_taxa <-
  ggplot() +
    # original had bubble dots - shape?
    {
      if (!is.null(.taxa)) {
        geom_point(data = taxa_amt_fin,
                   aes(x = x,
                       y = y,
                       color = taxa),
                   alpha = 0.5)
      } else {
        geom_point(data = taxa_amt_fin,
                   aes(x = x,
                       y = y),
                   alpha = 0.5,
                   colour = "red")
      }
    } +
    # scale_x_log10(limits = c(1, 100),
    #               labels = fancy_scientific,
    #               minor_breaks = xticks) +
    # scale_y_log10(limits = c(1e-5, 1),
    #               labels = fancy_scientific,
    #               minor_breaks = yticks) +
    labs(
      title = paste("Non Parametric Bootstrap:",
                    "Variation in Taxa coefficients"),
      x     = "Mean coefficient value",
      y     = "Coefficient of Variation"
    ) +
    theme_bw()
  
  plots <- list(plt_pig_r, plt_taxa)
  
  if (verbose) {
    print(plots)
  }
  
  # pigment ratio output
  pig_r_avg[indx]    <- pig_r_fin$x 
  pig_r_avg_sd[indx] <- pig_r_fin$std_dev
  
  # taxa contribution output
  taxa_amt_avg       <-
    matrix(taxa_amt_fin$x,  ncol = pig_r_row, nrow = df_row)
  taxa_amt_sd        <-
    matrix(taxa_amt_fin$std_dev,  ncol = pig_r_row, nrow = df_row)
  
  # ---- results ----
  results            <- list()
  results$replicates <- .nrep
  results$b          <- pig_r_avg
  results$b_sd       <- pig_r_avg_sd
  results$a          <- taxa_amt_avg
  results$a_sd       <- taxa_amt_sd
  results$plots      <- plots
  results$logs       <- logs
  results$ranges     <- list(pigment_range      = pig_rep_range, 
                             taxa_amount_ranges = taxa_amt_rep_range)
  
  return(results)
  # ---- end ----
}


