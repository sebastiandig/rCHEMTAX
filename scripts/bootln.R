bootln <- function(s, ssd, f, fsd)
################################################################################
#                                                                              # 
#               Plot the effect of log normal perturbation of s                #
#                                                                              #    
################################################################################
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
# Original: 2010-03-14  Matalb7  W.Whiten
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
root <- rprojroot::find_rstudio_root_file()
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
  ni   <- prod(dim(indx))
  tf   <- matrix(0,nrep,ni)
  tc   <- matrix(0,ns*nt)


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
  
  for (i in 1:nrep) {
    ss <- exp(
      matrix(rnorm(length(s1)),nrow(s1)) * sd1 + s1 # could be ncol, dont know the expected dimensions
    )
    # references another function
    # TODO: source(paste0(root,"/scripts/nnmatfactsd.R"))
    temp   <- nnmatfactsd(ss,ssd,f,fsd,info = data.frame(printitr = 1e6))
    cc     <- temp$cc
    ff     <- temp$ff
    info   <- temp$info

    tf[i,] <- Conj(ff[indx])
    tc[i,] <- Conj(cc)
    # not sure what toc is, but this is for printing only
    print(c(toc,info$rmsx,info$rmsxwt,info$itr/1000,info$conv*1e6))
  }

# disp(' ')
# disp('Range of variation in tf and tc')
# tf1=tf-repmat(mean(tf),nrep,1);
# disp([min(min(tf1)),max(max(tf1))])
# tc1=tc-repmat(mean(tc),nrep,1);
# disp([min(min(tc1)),max(max(tc1))])

  tf1 <- tf- pracma::repmat(mean(tf),nrep,1)
  tc1 <- tc- pracma::repmat(mean(tc),nrep,1)
  
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


# figure%
# %loglog(mean(tc)',(std(tc))','.')
# loglog(mean(tc)',(std(tc)./mean(tc))','.')
# title('Parametric Bootstrap, log normal:  Variation in c coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
# 
# return
# end
