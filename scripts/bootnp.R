bootnp <- function(s,ssd,f,fsd) {
################################################################################
#                                                                              # 
#           Plot the effect of non parametric bootstrap of s                   #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
#
# --------
# INPUTS:
# --------
# s       =  Matrix to be factored as a*b (best if a larger than b)
# ssd     = Matrix of standard deviations for x
# f       = Stabilising value for b, non zero locations, & initial value
# fsd     = Matrix of standard deviations for b
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
# Sebastian Di Geronimo (Fri Jun 10 12:01:31 2022)
  
  
# function bootnp(s,ssd,f,fsd)
# % bootnp  Plot the effect of non parametric bootstrap of s
# %  2010-02-27  Matlab7  W.Whiten
# %
# %  bootnp(s,sds,f,fsd)
# %  s    Matrix to be factored as a*b (best if a larger than b)
# %  ssd  Matrix of standard deviations for x
# %  f    Stabilising value for b, non zero locations, & initial value
# %  fsd  Matrix of standard deviations for b
#

  require("ggplot2")
  require("tictoc")
  require("pracma")
  require("rprojroot")
  
  root <- rprojroot::find_rstudio_root_file()
  nrep <-  10  
  # nrep=10;

  # references another function
  source(paste0(root,"/scripts/nnmatfactsd.R"))
  

  # [ns,np]=size(s); %#ok<NASGU>
  # [nt,np]=size(f);  %#ok<NASGU>
  # indx=find(f>0);
  # ni=length(indx);
  # tf=zeros(nrep,ni);
  # tc=zeros(nrep,ns*nt);

  # cc1=zeros(ns,nt);

  ns   <- dim(s)[1] # row of s
  np   <- dim(s)[2] # col of s
  nt   <- dim(f)[1] # row of f
  np   <- dim(f)[2] # col of f
  indx <- which(f > 0) 
  ni   <- length(indx)
  # ni   <- prod(dim(indx)[1])  # probably wont work
  tf   <- matrix(0, nrep, ni)
  tc   <- matrix(0, nrep, ns*nt)  

  cc1  <- matrix(0,ns,nt) 


  
# disp(' ')
# disp('   Time      rmsx      rmsxwt   itr/1000  conv*1e6')
# tic
# for i=1:nrep
# ind=ceil(rand(ns,1)*ns);
# [cc,ff,info]=nnmatfactsd(s(ind,:),ssd(ind,:),f,fsd,  ...
#                          struct('printitr',1e6));
# 
# tf(i,:)=ff(indx)';
#     cc1(:,:)=NaN;
#     cc1(ind,:)=cc;
#     tc(i,:)=cc1(:)';
# disp([toc,info.rmsx,info.rmsxwt,info.itr/1000,info.conv*1e6])
# end# disp(' ')
  
  tictoc::tic.clearlog()
  start <- tictoc::tic()
  
  for (i in 1:nrep) {
    ind <- pracma::ceil(pracma::rand(ns,1) * ns)
    temp   <- nnmatfactsd(s[ind,], ssd[ind,], f, fsd, 
                          info = list(printitr = 1e6))
    cc     <- temp$a
    ff     <- temp$b
    info   <- temp$info
    
    # tf[i,] <- Conj(ff[indx])
    tf[i,] <- t(ff[indx])
    cc1[,] <- NA
    cc1[ind,] <- cc
    # tc[i,] <- Conj(cc1)
    tc[i,] <- t(cc1)
    
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
  # disp([min(min(tc1(~isnan(tc1)))),max(max(tc1(~isnan(tc1))))])
  
  tf1 <- tf - pracma::repmat(apply(tf, 2, mean, na.rm = T),nrep,1)
  tc1 <- tc - pracma::repmat(apply(tc, 2, mean, na.rm = T),nrep,1)
  
  print('Range of variation in tf and tc')
  print(c(min(tf1, na.rm = T),max(tf1, na.rm = T)))
  print(c(min(tc1, na.rm = T),max(tc1, na.rm = T)))
  
  # figure
  # %loglog(mean(tf),std(tf),'o')
  # loglog(mean(tf),std(tf)./mean(tf),'o')
  # title('Non parametric bootstrap: Variation in f coefficients')
  # xlabel('Mean coefficient value')
  # ylabel('standard deviation/mean')
  
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
  df <- data.frame(x = apply(tf, 2, mean, na.rm=TRUE),
                   y = apply(tf, 2, function(x) 
                     sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE))
                   )
  
  print(ggplot() +
    # original had bubble dots - shape?
    geom_point(data = df,
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
      title = "Non parametric bootstrap: Variation in f coefficients",
      x     = "Mean coefficient value",
      y     = "standard deviation/mean"
    ) +
    theme_bw())

  # I don't think I need this chunk, but maybe
# mtc=zeros(ns*nt,1);
# stc=zeros(ns*nt,1);
# for i=1:ns*nt
#   t=tc(~isnan(tc(:,i)),i);
#   mtc(i)=mean(t);
#   stc(i)=std(t);
# end
  
  df2 <- data.frame(x = apply(tc, 2, mean, na.rm = TRUE),
                    y = apply(df, 2, function(x)
                      sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))
                    )

# figure
# %loglog(mtc',stc,'.')
# loglog(mtc,stc./mtc,'.')
# title('Non parametric bootstrap: Variation in c coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
# 
# return
# end

  print(ggplot() +
    # original had bubble dots - shape?
    geom_point(data = df2,
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
      title = "Non parametric bootstrap: Variation in c coefficients",
      x     = "Mean coefficient value",
      y     = "standard deviation/mean"
    ) +
    theme_bw())
}


