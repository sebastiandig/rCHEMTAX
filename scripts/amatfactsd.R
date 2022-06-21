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
  
  # brokewest -> chemtaxbrokewest() -> line 57 -> nnmatfactsd -> amatfactsd
  # amatfactsd(x,sdx,b) from nnmatfactsd
  # TODO: delete later 
  # x <- x
  # sdx <- sdx
  # b <- b
  
  ns <- dim(x)[1]
  nt <- dim(b)[1]
  
  a <- matrix(0, ns, nt)
  
  for (j in seq(ns)) {
  #  aa=lsqnonneg(b'./    repmat(sdx(j,:)',1,nt),
  #               (x(j,:)./sdx(j,:))');
    # sdx(j,:) is row j all columns, repeate
    # aa <- FUN( Conj(b)/ pracma::repmat(Conj(sdx[j,]),1,nt) , Conj(x[j,] / sdx[j,])  )
    
  # aa <- pracma::lsqnonneg(Conj(b)/ pracma::repmat(Conj(sdx[j,]),1,nt), 
  #                         Conj(x[j,] / sdx[j,])) 
    
  aa <- pracma::lsqnonneg(t(b)/ pracma::repmat(as.matrix((sdx[j,])),1,nt), 
                          as.vector(t((x[j,] / sdx[j,])))
                          )
  
  # a(j,:)=aa';
  a[j,] <- t(aa$x)
  }
  
  a
  
}