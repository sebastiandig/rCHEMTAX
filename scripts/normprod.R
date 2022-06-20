normprod <- function(s,c,f) {
################################################################################
#                                                                              # 
#               Normalise matrix product on last row of s & f                  #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
#
# --------
# INPUTS:
# --------
# s       = Original product matrix
# c       = Left factor of s
# f       = Right factor of s
#
# --------
# OUTPUTS:
# --------
# ss      = Matrix s rescaled for one in last column
# cc      = Left factor after scaling
# ff      = Right factor after scaling
# rms     = Root mean square of ss-cc*ff
#
# --------
# NOTES:
# --------
# Original: 2010-03-28  Matlab7  W.Whiten
#
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sun Jun 19 14:53:10 2022)
 
  # function [ss,cc,ff,rms]=normprod(s,c,f)
  # % normprod  Normalise matrix product on last row of s & f
  # %  2010-03-28  Matlab7  W.Whiten
  # %
  # % [ss,cc,ff,rms]=normprod(s,c,f)
  # %  s   Original product matrix
  # %  c   Left factor of s
  # %  f   Right factor of s
  # %
  # %  ss  Matrix s rescaled for one in last column
  # %  cc  Left factor after scaling
  # %  ff  Right factor after scaling
  # %  rms Root mean square of ss-cc*ff
  
  library("pracma")
  
  # [ns,np]=size(s);
  # nt=size(f,1);
  
  ns <- dim(s)[1]
  np <- dim(s)[2]
  nt <- dim(f)[1]
  
  # s1=1./(s(:,end)+1e-100);
  # ss=s.*repmat(s1,1,np);
  
  s1 <- as.matrix(1 / (s[,ncol(s)] + 1e-100))
  ss <- s * pracma::repmat(s1, 1, np)
  
  # s2=f(:,end);
  # ff=f.*repmat(1./(s2+1e-100),1,np);
  
  s2 <- as.matrix(f[,ncol(f)])
  ff <- f * pracma::repmat(1/(s2 + 1e-100), 1, np)
  
  # t=repmat(s1,1,nt).*c;
  # cc=t.*repmat(s2',ns,1);

  t  <- pracma::repmat(s1, 1, nt) * c
  cc <- t * pracma::repmat(t(s2), ns, 1)
  
  # % get rms
  # t=ss-cc*ff;
  # rms=sqrt(sum(t(:).^2)/length(t(:)));
  
  t   <- ss - cc %*% ff 
  # rms <- sqrt( sum(t^2)  / length(t))
  rms <- sqrt(mean(t^2))
  
  # return
  # end 
  
  results <- list(ss = ss,
                  cc = cc,
                  ff = ff,
                  rms = rms)
  
  return(results)
}
