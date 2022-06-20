amatfactsd <- function(x, sdx, b){
################################################################################
#                                                                              # 
#               Calculate Factor a in min(x-ab) using least square             #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
#
# --------
# INPUTS:
# --------
# x       = Matrix to be factored as a*b
# sdx     = Standard deviations of x
# b       = Value of b
#
# --------
# OUTPUTS:
# --------
# a       = Result left factor of x
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
# Sebastian Di Geronimo (Thu Jun 02 18:30:06 2022)
  

  
# amatfact  Calculate factor a in min(x-ab) using lsqnonneg
# 
# ns=size(x,1);
# nt=size(b,1);
# 
# a=zeros(ns,nt);
# for j=1:ns
# %        aa=lsqnonneg(bb,x(j,:)'); 
#     aa=lsqnonneg(b'./repmat(sdx(j,:)',1,nt),(x(j,:)./sdx(j,:))');
#     a(j,:)=aa';
# end
# 
# return
# end

  library("pracma")
  x <- matrix(c(0.0372,    0.2869,
                0.6861,   0.7071,
                0.6233,    0.6245,
                0.6344,    0.6170),
              nrow = 4,
              ncol = 2)
  
  b <- c(0.8587,
              0.1781,
              0.0747,
              0.8405)
  
  ns <- dim(x)[1]
  nt <- length(b)
  
  a <- matrix(0, ns, nt)
  
  for (j in seq(ns)) {
    
  # TODO: figure out the dimensions of sdx ----------------------------------------
    
  #  aa=lsqnonneg(b'./    repmat(sdx(j,:)',1,nt),
  #               (x(j,:)./sdx(j,:))');
    # sdx(j,:) is row j all columns, repeate
    # aa <- FUN( Conj(b)/ pracma::repmat(Conj(sdx[j,]),1,nt) , Conj(x[j,] / sdx[j,])  )
    
  # aa <- pracma::lsqnonneg(Conj(b)/ pracma::repmat(Conj(sdx[j,]),1,nt), 
  #                         Conj(x[j,] / sdx[j,]))  
  aa <- pracma::lsqnonneg(t(b)/ pracma::repmat(t(sdx[j,]),1,nt), 
                          t(x[j,] / sdx[j,]))
  
  # a(j,:)=aa';
  a[j,] <- t(aa)
  }
  
}