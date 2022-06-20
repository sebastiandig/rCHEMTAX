randstart <- function(s,ssd,f,fsd) {
################################################################################
#                                                                              # 
#   Plot the effect of random starts for s% effect of random start locations   #
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
# NA
#
# --------
# NOTES:
# --------
# Original: 2010-02-26  Matlab7  W.Whiten
#  
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sun Jun 19 14:01:35 2022)  

  
  # function randstart(s,ssd,f,fsd)
  # % randstart  Plot the effect of random starts for s% effect of random start locations
  # %  2010-02-26  Matlab7  W.Whiten
  # %
  # %  regplot(s,sds,f,fsd)
  # %  s    Matrix to be factored as a*b (best if a larger than b)
  # %  ssd  Matrix of standard deviations for x
  # %  f    Stabilising value for b, non zero locations, & initial value
  # %  fsd  Matrix of standard deviations for b
  
  library("tictoc")
  library("pracma")
  library("ggplot2")
  
  # set directory for saving 
  root <- rprojroot::find_rstudio_root_file()
  
  # references another function
  source(paste0(root,"/scripts/nnmatfactsd.R"))
  
  # nrep=10;
  nrep <- 10 
  
  # [ns,np]=size(s); %#ok<NASGU>
  # [nt,np]=size(f); %#ok<NASGU>
  # indx=find(f>0);
  # ni=length(indx);
  # tf=zeros(nrep,ni);
  # tc=zeros(nrep,ns*nt);
  
  ns   <- dim(s)[1] # row of s
  np   <- dim(s)[2] # col of s
  nt   <- dim(f)[1] # row of f
  np   <- dim(f)[2] # col of f
  indx <- which(f > 0)
  ni   <- length(indx)
  tf   <- matrix(0, nrep, ni)
  tc   <- matrix(0, nrep, ns*nt)
  
  
  # disp('   Time      rmsx      rmsxwt   itr/1000  conv*1e6')
  # tic
  # for i=1:nrep
  # [cc,ff,info]=nnmatfactsd(s,ssd,f,fsd,  ...
  #                          struct('inita',rand(ns,nt),'maxitr',30000,'printitr',1e12,   ...
  #                                 'conv',1e-10));
  #
  # tf(i,:)=ff(indx)';
  # tc(i,:)=cc(:)';
  # disp([toc,info.rmsx,info.rmsxwt,info.itr/1000,info.conv*1e6])
  # end
  
  tictoc::tic.clearlog()
  start <- tictoc::tic()
  
  for (i in 1:nrep) {
    temp   <- nnmatfactsd(s,
                          ssd,
                          f,
                          fsd,
                          info = list(
                            inita = pracma::rand(ns, nt),
                            maxitr = 30000,
                            printitr = 1e12,
                            conv = 1e-10
                          ))
    cc     <- temp$cc
    ff     <- temp$ff
    info   <- temp$info
    
    tf[i,] <- t(ff[indx])
    tc[i,] <- t(cc)
    
    end <- (tictoc::toc(quiet = T))$toc - start
    print(c(Time   = end,
            rmsx   = info$rmsx,
            rmsxwt = info$rmsxwt,
            itr    = info$itr/1000,
            conv   = info$conv*1e6))
  }
  
  # disp(' ')
  # disp('Range of variation in tf and tc')
  # tf1=tf-repmat(mean(tf),10,1);
  # disp([min(min(tf1)),max(max(tf1))])
  # tc1=tc-repmat(mean(tc),10,1);
  # disp([min(min(tc1)),max(max(tc1))])
  
  # may need to wrap mean(tf/tc) -> as.matrix(apply(tf, 1/2, mean))
  tf1 <- tf - pracma::repmat(as.matrix(mean(tf, na.rm = T)),nrep,1)
  tc1 <- tc - pracma::repmat(as.matrix(mean(tc, na.rm = T)),nrep,1)
  
  print('Range of variation in tf and tc')
  print(c(min(tf1, na.rm = T),max(tf1, na.rm = T)))
  print(c(min(tc1, na.rm = T),max(tc1, na.rm = T)))
  
  
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
  
  
  # figure
  # loglog(mean(tf),std(tf)./mean(tf),'o')
  # title('Variation in f coefficients')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
  # tf
  df <- data.frame(x = apply(tf, 2, mean, na.rm=TRUE),
                   y = apply(tf, 2, function(x) 
                     sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
  )
  
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
    labs(title = "Variation in f coefficients",
         x     = "Mean coefficient value",
         y     = "standard deviation/mean") +
    theme_bw()
  
  
  # figure
  # loglog(mean(tc)',(std(tc)./mean(tc))','.')
  # title('Variation in c coefficients')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
  # tc
  # might have issue with transpose?
  df2 <- data.frame(x = apply(tc, 2, mean, na.rm = TRUE),
                    y = apply(tc, 2, function(x)
                      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
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
    labs(title = "Variation in c coefficients",
         x     = "Mean coefficient value",
         y     = "standard deviation/mean") +
    theme_bw()
  
}