bootln <- function(s, ssd, f, fsd)
################################################################################
#                                                                              # 
#               Plot the effect of log normal perturbation of s                #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
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
# Currently none as only plots
#
# --------
# NOTES:
# --------
# Original: 2010-02-27  Matlab7  W.Whiten
#
# --------
# References:
# --------
#
# --------
# Author:
# --------
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
# 

  require("ggplot2")
  require("tictoc")
  require("pracma")

  root <- rprojroot::find_rstudio_root_file()
  nrep <-  10
# nrep=10;
# 
# [ns,np]=size(s); %#ok<NASGU>
# [nt,np]=size(f); %#ok<NASGU>
# indx=find(f>0);
# ni=length(indx);
# tf=zeros(nrep,ni);
# tc=zeros(nrep,ns*nt);
# 

  ns   <-  dim(s)[1] # row of s
  np   <-  dim(s)[2] # col of s
  nt   <-  dim(f)[1] # row of f
  np   <-  dim(f)[2] # col of f
  indx <- which(f > 0)
  ni   <- length(indx)
  # ni   <- prod(dim(indx))
  tf   <- matrix(0, nrep, ni)
  tc   <- matrix(0, nrep, ns*nt)


# % parameters for log normal distribution
# s1=max(s,ssd/2);
# sd1=log(1+(ssd./s1).^2);
# s1=log(s1)-sd1/2;
# sd1=sqrt(sd1);
# 

  s1   <- pmax(s, ssd/2)
  sd1  <- log(1+(ssd/s1)^2)
  s1   <- log(s1)-sd1/2
  sd1  <- sqrt(sd1)
  
# disp(' ')
# disp('   Time      rmsx      rmsxwt   itr/1000  conv*1e6')
# tic
# for i=1:nrep
# ss=exp(randn(size(s1)).*sd1+s1);
# [cc,ff,info]=nnmatfactsd(ss,ssd,f,fsd,struct('printitr',1e6));
# 
# tf(i,:)=ff(indx)';
#     tc(i,:)=cc(:)';
# disp([toc,info.rmsx,info.rmsxwt,info.itr/1000,info.conv*1e6])
# end
# 
  tictoc::tic.clearlog()
  start <- tictoc::tic()
  
  # references another function nnmatfactsd
  # TODO: source(paste0(root,"/scripts/nnmatfactsd.R"))
  
  for (i in 1:nrep) {
    tictoc::tic()
    ss <- exp(
      # runif over rnorm because rnorm makes negative numbers
      matrix(runif(length(s1)),nrow(s1)) * sd1 + s1 # could be ncol, dont know the expected dimensions
    )

    temp   <- nnmatfactsd(ss, ssd, f, fsd, info = data.frame(printitr = 1e6))
    cc     <- temp$cc
    ff     <- temp$ff
    info   <- temp$info

    tf[i,] <- Conj(ff[indx])
    tc[i,] <- Conj(cc)

    end <- (tictoc::toc(quiet = T))$toc - start
    print(c(Time   = end,
            rmsx   = info$rmsx,
            rmsxwt = info$rmsxwt,
            itr    = info$itr/1000,
            conv   = info$conv*1e6))
  }

# disp(' ')
# disp('Range of variation in tf and tc')
# tf1=tf-repmat(mean(tf),nrep,1);
# disp([min(min(tf1)),max(max(tf1))])
# tc1=tc-repmat(mean(tc),nrep,1);
# disp([min(min(tc1)),max(max(tc1))])

  tf1 <- tf - pracma::repmat(mean(tf),nrep,1)
  tc1 <- tc - pracma::repmat(mean(tc),nrep,1)
  
  print('Range of variation in tf and tc')
  print(c(min(tf1),max(tf1)))
  print(c(min(tc1),max(tc1)))
  
# figure
# %loglog(mean(tf),std(tf),'o')
# loglog(mean(tf),std(tf)./mean(tf),'o')
# title('Parametric Bootstrap, log normal:  Variation in f coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
#
  
  # allow loglog style plotting from Matlab
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
  
  # tf
  df <- data.frame(x = apply(tf, 2, mean),
                   y = apply(tf, 2, function(x) sd(x)/mean(x)))
  
  ggplot() +
    # original had bubble dots - shape?
    geom_point(data = df,
               aes(x = x,
                   y = y),
               colour = "red") +
    scale_x_log10(limits = c(1, 100),
                  labels = fancy_scientific,
                  minor_breaks = xticks) +
    scale_y_log10(limits = c(1e-5, 1),
                  labels = fancy_scientific,
                  minor_breaks = yticks) +
    labs(title = "Parametric Bootstrap, log normal:  Variation in f coefficients",
         x     = "Mean coefficient value",
         y     = "standard deviation/mean") +
    theme_bw()
  
# figure%
# %loglog(mean(tc)',(std(tc))','.')
# loglog(mean(tc)',(std(tc)./mean(tc))','.')
# title('Parametric Bootstrap, log normal:  Variation in c coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
# 
# return
# end

  # tc
  # might have issue with transpose?
  df2 <- data.frame(x = apply(tc, 2, mean),
                   y = apply(tc, 2, function(x) sd(x)/mean(x))
                   )
  
  ggplot() +
    # original had bubble dots - shape?
    geom_point(data = df2,
               aes(x = x,
                   y = y),
               colour = "red") +
    scale_x_log10(limits = c(1, 100),
                  labels = fancy_scientific,
                  minor_breaks = xticks) +
    scale_y_log10(limits = c(1e-5, 1),
                  labels = fancy_scientific,
                  minor_breaks = yticks) +
    labs(
      title = "Parametric Bootstrap, log normal:  Variation in c coefficients",
      x = "Mean coefficient value",
      y = "standard deviation/mean"
    ) +
    theme_bw()