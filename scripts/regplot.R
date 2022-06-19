regplot <- function(s,ssd,f,fsd) {
################################################################################
#                                                                              # 
#          Plot the effect of weighting of f in positive matrix factors        #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# Ploting function
#
# --------
# INPUTS:
# --------
# s       = Matrix to be factored as a*b (best if a larger than b)
# ssd     = Matrix of standard deviations for x
# f       = Stabilising value for b, non zero locations, & initial value
# fsd     = Matrix of standard deviations for b
#
# --------
# OUTPUTS:
# --------
# Plot
#
# --------
# NOTES:
# --------
# Original: 2010-02-21  Matlab7  W.Whiten
#
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sat Jun 18 20:54:16 2022)

  # function regplot(s,ssd,f,fsd)
  # % regplot  Plot the effect of weighting of f in positive matrix factors
  # %  2010-02-21  Matlab7  W.Whiten
  # %
  # %  regplot(s,sds,f,fsd)
  # %  s    Matrix to be factored as a*b (best if a larger than b)
  # %  ssd  Matrix of standard deviations for x
  # %  f    Stabilising value for b, non zero locations, & initial value
  # %  fsd  Matrix of standard deviations for b
  
  # % plot the effect of prior on f
  # t=[10, sqrt(10)];
  # fsdx=[ t 0.1*t 0.01*t 0.02 0.001*t 0.0001*t]';
  # n=length(fsdx);
  # rmsx=zeros(n,1);
  # rmsxwt=zeros(n,1);
  
  library("tictoc")
  library("ggplot2")
  
  # set directory for saving 
  root <- rprojroot::find_rstudio_root_file()
  
  # source(paste0(root, "scripts/nnmatfactsd.R"))
  
  t      <- c(10, sqrt(10))
  fsdx   <- c(t, 0.1*t, 0.01*t, 0.02, 0.001*t, 0.0001*t)
  n      <- length(fsdx)
  rmsx   <- matrix(0, n, 1)
  rmsxwt <- matrix(0, n, 1)
  
  # for i=1:n
  #     # disp(' ')
  #     disp(['Factor=' num2str(fsdx(i))])
  #     tic;[cc,ff,info]=nnmatfactsd(s,ssd,f,fsdx(i)*fsd);toc %#ok<ASGLU>
  #     rmsx(i)=info.rmsx;
  #     rmsxwt(i)=info.rmsxwt;
  # end
  
  for (1 in 1:n) {
    
    print(paste("Factor = ", fsdx[i]))
    
    tictoc::tic()
    temp      <- nnmatfactsd(s,ssd,f,fsdx[i]*fsd)
    
    info      <- temp$info
    rmsx[i]   <- info$rmsx
    rmsxwt[i] <- info$rmsxwt
    
    tictoc::toc()
  }    
  
  # figure
  # semilogx(fsdx,rmsxwt)
  # title('Weighted root mean square error vs sd factor')
  # xlabel('Factor on f sd')
  # ylabel('Root mean square')
  
  # allow y log style plotting from Matlab in semilogy
  fancy_scientific <- function(l) {
    # turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    # quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    # turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    # return this as an expression
    parse(text=l)
  }
  
  yticks = outer((1:10),(10^(-5:-1)))
  xticks = outer((1:10),(10^(0:1)))
  
  
  df <- data.frame(fsdx = fsdx,rmsxwt= rmsxwt)
  
  ggplot(df, aes(x = fsdx, y = rmsxwt)) +
    geom_path() +
    labs(title = 'Weighted root mean square error vs sd factor',
         x = 'Factor on f sd',
         y ='Root mean square') +
    scale_x_log10(limits = c(1, 100),
                  labels = fancy_scientific,
                  minor_breaks = xticks) +
    theme_bw()
  
  # return
  # end
  
  
}

