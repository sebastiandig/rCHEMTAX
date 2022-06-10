bootnp <- function(s,ssd,f,fsd) {
################################################################################
#                                                                              # 
#           Plot the effect of non parametric bootstrap of s                   #
#                                                                              #    
################################################################################
# INPUTS:
# --------
# s       =  Matrix to be factored as a*b (best if a larger than b)
# ssd = Matrix of standard deviations for x
# f = Stabilising value for b, non zero locations, & initial value
# fsd = Matrix of standard deviations for b
# --------
# OUTPUTS:
# --------
#  =
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
# Author (Fri Jun 10 12:01:31 2022)
  
  
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
# nrep=10;
# %nrep=20;
# 
# [ns,np]=size(s); %#ok<NASGU>
# [nt,np]=size(f);  %#ok<NASGU>
# indx=find(f>0);
# ni=length(indx);
# tf=zeros(nrep,ni);
# tc=zeros(nrep,ns*nt);
# 
# cc1=zeros(ns,nt);
# 
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
# end
# 
# disp(' ')
# disp('Range of variation in tf and tc')
# tf1=tf-repmat(mean(tf),nrep,1);
# disp([min(min(tf1)),max(max(tf1))])
# tc1=tc-repmat(mean(tc),nrep,1);
# disp([min(min(tc1(~isnan(tc1)))),max(max(tc1(~isnan(tc1))))])
# 
# figure
# %loglog(mean(tf),std(tf),'o')
# loglog(mean(tf),std(tf)./mean(tf),'o')
# title('Non parametric bootstrap: Variation in f coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
# 
# mtc=zeros(ns*nt,1);
# stc=zeros(ns*nt,1);
# for i=1:ns*nt
# t=tc(~isnan(tc(:,i)),i);
# mtc(i)=mean(t);
# stc(i)=std(t);
# end
# figure
# %loglog(mtc',stc,'.')
# loglog(mtc,stc./mtc,'.')
# title('Non parametric bootstrap: Variation in c coefficients')
# xlabel('Mean coefficient value')
# ylabel('standard deviation/mean')
# 
# return
# end
}