normprod <- function(.df,c,f) {
################################################################################
#                                                                              # 
#               Normalize matrix product on last row of s & f                  #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: -------
#
# ---- INPUTS: -----------
# df = Original product matrix
# c = Left factor of df
# f = Right factor of df
#
# ---- OUTPUTS: ----------
# ss = Matrix .df rescaled for one in last column
# cc = Left factor after scaling
# ff = Right factor after scaling
# rms = Root mean square of ss-cc*ff
#
# ---- NOTES: ------------
# Original: 2010-03-28  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sun Jun 19 14:53:10 2022)
 
  # %  2010-03-28  Matlab7  W.Whiten
  
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
  
  ns <- dim(.df)[1]
  np <- dim(.df)[2]
  nt <- dim(f)[1]
  
  # s1=1./(s(:,end)+1e-100);
  # ss=s.*repmat(s1,1,np);
  
  s1 <- as.matrix(1 / (.df[,ncol(.df)] + 1e-100))
  ss <- .df * pracma::repmat(s1, 1, np)
  
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
