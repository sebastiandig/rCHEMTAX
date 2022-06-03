bmatfactsd <- function(x, sdx, a, b0, sdb) {
  
################################################################################
#                                                                              # 
#       Calculate factor b in min(x-ab, b-b0) using Least Linear Squares       #
#                                                                              #    
################################################################################
# INPUTS:
#--------
# input = 
#        x    Matrix to be factored as a*b
#        sdx  Standard deviations of x
#        a    result of a=amatfactsd(x,sdx,b,sdb)
#        b0   Initial estimate for b
#        sdb  Standard deviations of b0?
#        
##########
# OUTPUTS:
#---------
# output =  
#         b    Result left factor of x
##########
# -----Notes:-----
# Original: 2010-03-14  Matalb7  W.Whiten
# -----References:-----
#
# -----Author:-----
# Sebastian Di Geonimo (Thu Jun 02 19:54:42 2022)
  

# function b=bmatfactsd(x,sdx,a,b0,sdb)
# % bmatfact  Calculate factor b in min(x-ab,b-b0) using lsqnonneg
# %  
# %
# % a=amatfactsd(x,sdx,b,sdb)
# %  x    Matrix to be factored as a*b
# %  sdx  Standard deviations of x
# %  b    Value of b
# %  b0   Initial estimate for b
# %  
# %
# %  a    Result left factor of x
# 
# %ns=size(x,1);
# [nt,np]=size(b0);
# 
# indx1=b0~=0;
# b=zeros(nt,np);
# for j=1:np
# indx=indx1(:,j);
# %        t=lsqnonneg(a(:,indx),x(:,j));
# t=lsqnonneg([a(:,indx)./  ...
#              repmat(sdx(:,j),1,sum(indx));diag(1./sdb(indx,j))], ...
#             [x(:,j)./sdx(:,j);b0(indx,j)./sdb(indx,j)]);
# b(indx,j)=t;
# end
# 
# return
# end
  
  library("pracma")
  
  # dimentions of b0
  b0_dim <- dim(b0)
  nt <- b0_dim[1] # row
  np <- bo_dim[2] # col
  
  # index where b0 does not equal 0
  indx1 <- b0 != 0
  
  # create empty matrix of 0s with nt and np dimenstions
  b <- matrix(0, nt, np)
  
  for (j in seq(np)) {
    indx <- indx1[,j]
    
    # t=lsqnonneg([a(:,indx)./ repmat(sdx(:,j),1,sum(indx));
    #              diag(1./sdb(indx,j))], 
    #             
    #              [x(:,j)./sdx(:,j);
    #               b0(indx,j)./sdb(indx,j)])
    
    t <- pracma::lsqnonneg(
      
      cbind( # might be c or rbind, not sure of input dimensions
        a[,indx] / pracma::repmat(sdx[,j],1, sum(indx)),
        diag(1 / sdb[indx,j])
      ),
      cbind( # might be c or rbind, not sure of input dimensions
        (x[,j] / sdx[,j]), 
        (b0[indx,j] / sdb[indx,j])
      )
    )            
  }
}