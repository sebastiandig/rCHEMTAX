bootln <- function(.df, .df_sd, .pig_r, .pig_r_sd, .info = NULL, .nrep = 10, 
                   verbose = TRUE, .pigm = NULL, .taxa = NULL){
################################################################################
#                                                                              # 
#               Plot the effect of log normal perturbation of s                #
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
# Currently none as only plots
#
# ---- NOTES: ------------
# Original: 2010-02-27  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Thu Jun 09 16:18:57 2022)
 
  # function bootln(s,ssd,f,fsd)
  # % bootln  Plot the effect of log normal perturbation of s
  # %  2010-02-27  Matlab7  W.Whiten
  # %
  # %  bootln(s,sds,f,fsd)
  # %  s    Matrix to be factored as a*b (best if a larger than b)
  # %  ssd  Matrix of standard deviations for x
  # %  f    Stabilising value for b, non zero locations, & initial value
  # %  fsd  Matrix of standard deviations for b

  
  # ---- TODO: remove after testing ----
  # .df       = s
  # .df_sd    = ssd
  # .pig_r    = f0
  # .pig_r_sd = fsd
  # .nrep = 5
  # .info = NULL
  # verbose = T
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
  # ni         <- prod(dim(indx))
  pig_r_rep    <- matrix(0, .nrep, n_pig)
  taxa_amt_rep <- matrix(0, .nrep, df_row*pig_r_row)

  # initialize empty matrix
  pig_r_avg    <- matrix(0, pig_r_row, pig_r_col)
  pig_r_avg_sd <- matrix(0, pig_r_row, pig_r_col)
  taxa_amt_avg <- matrix(0, df_row, df_col)
  taxa_amt_sd  <- matrix(0, df_row, df_col)
  
  # ---- parameters for log normal distribution ----
  
  # s1=max(s,ssd/2);
  # sd1=log(1+(ssd./s1).^2);
  # s1=log(s1)-sd1/2;
  # sd1=sqrt(sd1);
  df_mx     <- pmax(.df, .df_sd / 2)
  df_sd     <- log(1 + (.df_sd / df_mx)^2)
  df_log    <- log(df_mx) - df_sd / 2
  df_sd_log <- sqrt(df_sd)
  
  # ---- options ----
  # create open list if no defaults are given for `info`  
  if (is.null(.info)) {
    .info <- list()
  }
  
  # set defaults options
  deflt <- list(
    printitr = 1e6
  )
  
  # add/replace default options if not set
  info  <- initstruct(.info, deflt)
  
  maxitr <- info$maxitr
  printitr <- info$printitr
  # convitr <- info$conitr
  
  # ---- initialize list for logs ----
  logs <- list() 
  

  # ss=exp(randn(size(s1)).*sd1+s1);
  # [cc,ff,info]=nnmatfactsd(ss,ssd,f,fsd,struct('printitr',1e6));
  # tf(i,:)=ff(indx)';
  # tc(i,:)=cc(:)';
  # disp([toc,info.rmsx,info.rmsxwt,info.itr/1000,info.conv*1e6])

  
  # ---- start loop ----
  tictoc::tic.clearlog()
  start <- tictoc::tic()

  for (i in seq(.nrep)) {
    tictoc::tic()
    
    cat(sprintf('\nRandom Start Number: %02d of %02d\n', i, .nrep))
    
    # runif over rnorm because rnorm makes negative numbers
    # matrix(runif(length(s1)),nrow(s1)) * sd1 + s1 # could be ncol, dont know 
    # the expected dimensions
    ss <- exp(
      as.matrix(pracma::rand(df_row, pig_r_col)) * df_sd_log + df_log # better implementation
    )

    temp          <- nnmatfactsd(ss, .df_sd, .pig_r, .pig_r_sd, 
                                 .info = info)
    taxa_amt_temp <- temp$a
    pig_r_temp    <- temp$b
    info          <- temp$info
    
    # log taxa contribution and pigment ratio out
    pig_r_rep[i,]    <- t(pig_r_temp[indx])
    taxa_amt_rep[i,] <- t(taxa_amt_temp)
    
    # log each run
    rep_num          <- paste0("rep_num_", i)
    logs[[rep_num]]  <- temp

    end <- (tictoc::toc(quiet = T))$toc - start

    # when printitr is above maxitr, will only print last iterations
    if (!(info$printitr <= info$maxitr) & verbose) {
      # cat(sprintf('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:     dRMS:   dTaxa RMS:   dPig RMS:\n'))
      cat('\nIter:    RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:    ',
          'dRMS:   dTaxa RMS:   dPig RMS:\n')
      cat(sprintf('%5i%#8.3g%#10.3g%#21.3g%#19.3g%#10.2e%13.2e%#12.2e\n',
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

  # disp(' ')
  # disp('Range of variation in tf and tc')
  # tf1=tf-repmat(mean(tf),nrep,1);
  # disp([min(min(tf1)),max(max(tf1))])
  # tc1=tc-repmat(mean(tc),nrep,1);
  # disp([min(min(tc1)),max(max(tc1))])
  
  # ---- display range of values for each pigment and taxa contributions ----
  pig_rep_range    <- 
    pig_r_rep - pracma::repmat(apply(pig_r_rep, 2, mean, na.rm = T),.nrep,1)
  df_pig_rep_range <- 
    taxa_amt_rep - pracma::repmat(apply(taxa_amt_rep, 2, mean, na.rm = T),.nrep,1)
  if (verbose) {
    
    cat('\nRange of variation in Pigment Ratios and Taxa Contribution:
      \n-----------------------------\n',
      sprintf('\nPigment Ratio:\n%14s%7.3f\n%14s%7.3f\n',
              "Max:", min(pig_rep_range, na.rm = T), 
              "Min:", max(pig_rep_range, na.rm = T)),
      sprintf('\nTaxa Contributon:\n%14s%7.3f\n%14s%7.3f\n',
              "Max:", min(df_pig_rep_range, na.rm = T), 
              "Min:", max(df_pig_rep_range, na.rm = T)),
      '\n-----------------------------\n')

  }
  
  # ---- calc avg and sd for final matrices ----
  # pig_r_rep
  pig_r_fin <-
    data.frame(
      x       = apply(pig_r_rep, 2, mean, na.rm = TRUE),
      std_dev = apply(pig_r_rep, 2, sd, na.rm = TRUE),
      y       = apply(pig_r_rep, 2, function(x)
                sd(x , na.rm = TRUE) / mean(x, na.rm = TRUE))
    )
  # TODO: work on getting names on graph
  # .pigm = pigm
  # 
  # if (!is.null(.taxa)) {
  #   pig_r_fin <- cbind(pig_r_fin,
  #                         pigm = rep(.pigm, pig_r_col)[indx])
  # }
  
  
  # taxa contribution
  taxa_amt_fin <- 
    data.frame(
      x       = apply(taxa_amt_rep, 2, mean, na.rm = TRUE),
      std_dev = apply(taxa_amt_rep, 2, sd, na.rm = TRUE),
      y       = apply(taxa_amt_rep, 2, function(x)
                      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
  )
  
  if (!is.null(.taxa)) {
    taxa_amt_fin <- cbind(taxa_amt_fin,
                          taxa = rep(.taxa, df_row))
  }
  
  

  
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
    #               minor_breaks = yticks_minor) +
    labs(title = paste("Parametric Bootstrap, log normal:",
                       "Variation in Pigment Ratio coefficients"),
         x     = "Mean coefficient value",
         y     = "Coefficient of Variation") +
    theme_bw()
  

  # taxa_amt_rep
  # plt_taxa <- 
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
    #               minor_breaks = yticks_minor) +
    labs(
      title = paste("Parametric Bootstrap, log normal:",
                    "Variation in Taxa coefficients"),
      x = "Mean coefficient value",
      y = "Coefficient of Variation"
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
  
  return(results)
}