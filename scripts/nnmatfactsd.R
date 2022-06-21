nnmatfactsd <- function(x,sdx,b0,sdb,info=NULL){
################################################################################
#                                                                              # 
#            Non negative matrix factors x=a*b with prior b0 & sd              #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# Can normalise matrices after function exits  eg  a*inv(d), d*b
#
# --------
# INPUTS:
# --------
# x       = Matrix to be factored as a*b (best if a larger than b)
# sdx     = Matrix of standard deviations for x
# b0      = Stabilising value for b, non zero locations, & initial value
# sdb     = Matrix of standard deviations for b
# info    = Information for calculation (to override default values)
#           .maxitr   Maximum number of iterations
#           .convitr  Iterations between convergence tests
#           .printitr Count of iterations between printing
#           .conv     Converge test value for rms change
#           .inita    Initial a value
#           .initb    Initial b value
#           .maxaitr  Maximum for initial iterations for a
#
# --------
# OUTPUTS:
# --------
# a       = Left factor of x
# b       = Right factor of x
# info    = Initial info and
#           .itr    Iterations used
#           .conv   Converge level for error
#           .rms    Total rms error weighted for x & b
#           .rmsx   Rms error in x unweighted
#           .rmdxwt Weighted rms error in x
#           .rmsb   Rms error in b unweighted
#           .rmsbwt Weighted rms error in b
#
# --------
# NOTES:
# --------
# Original: 2010-03-12  Matlab7  W.Whiten
#
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sun Jun 19 22:22:08 2022) 
  
  # % nnmatfactsd  Non negative matrix factors x=a*b with prior b0 & sd
  # %              extends algorithm of Lee & Seung
  # %  2010-03-12  Matlab7  W.Whiten
  # %
  # % [a,b,info]=nnmatfactsd(x,sdx,b0,sdb,info)
  # %  x    Matrix to be factored as a*b (best if a larger than b)
  # %  sdx  Matrix of standard deviations for x
  # %  b0   Stabilising value for b, non zero locations, & initial value
  # %  sdb  Matrix of standard deviations for b
  # %  info Information for calculation (to override default values)
  # %        .maxitr   Maximum number of iterations
  # %        .convitr  Iterations between convergence tests
  # %        .printitr Count of iterations between printing
  # %        .conv     Converge test value for rms change
  # %        .inita    Initial a value
  # %        .initb    Initial b value
  # %        .maxaitr  Maximum for initial iterations for a
  # %
  # %  a    Left factor of x
  # %  b    Right factor of x
  # %  info Initial info and
  # %        .itr    Iterations used
  # %        .conv   Converge level for error
  # %        .rms    Total rms error weighted for x & b
  # %        .rmsx   Rms error in x unweighted
  # %        .rmdxwt Weighted rms error in x
  # %        .rmsb   Rms error in b unweighted
  # %        .rmsbwt Weighted rms error in b
  # %
  # % Can normalise matrices after function exits  eg  a*inv(d), d*b
  # % 
  
  library("pracma")
  
  root <- rprojroot::find_rstudio_root_file()
  
  source(paste0(root,"/scripts/initstruct.R"))
  source(paste0(root,"/scripts/amatfactsd.R"))
  
  # remove later to info <- NULL
  # Inputs here x,sdx,b0,sdb,info=NULL
  # Input from brokewest.R -> nnmatfactsd(s,ssd,f0,fsd)
  # x <- s
  # sdx <- ssd
  # b0 <- f0
  # sdb <- fsd
  # info <- NULL
  
  # % default for info
  # if(nargin<5)
  #   info=struct();
  # end
  
  if (is.null(info)){
    info <- list()
  }
  
  # % matrix sizes
  # [ns,np]=size(x);
  # [nt,nx]=size(b0);
  # if(np~=nx); error('matfactsd: matrix size error'); end
  # nx=ns*np;
  # na=ns*nt;
  # nb0=sum(sum(b0~=0));
  # nerr=nx+nb0;
  
  ns   <- dim(x)[1]
  np   <- dim(x)[2]
  nt   <- dim(b0)[1]
  nx   <- dim(b0)[2]
  if (np != nx) stop('matfactsd: matrix size error\nnot the same size')
  
  nx   <- ns * np
  na   <- ns * nt
  nb0  <- sum(b0!=0)
  nerr <- nx + nb0
  
  # % transformed values as used in calculation
  # sdx=sdx+1e-100;
  # sdb=sdb+1e-100;
  # w2x=1./(sdx.^2);
  # w2b=1./(sdb.^2);
  # xw=x.*w2x;
  # b0w=b0.*w2b;
  
  sdx <- sdx + 1e-100
  sdb <- sdb + 1e-100
  w2x <- 1 / (sdx^2)
  w2b <- 1 / (sdb^2)
  xw  <- x * w2x
  b0w <- b0 * w2b
  
  # % replace default values by given values
  # deflt=struct('maxitr',20000,'convitr',500,'printitr',1000,  ...
  #              'conv',1e-6,'initb',b0,'maxaitr',-1);
  # info=initstruct(info,deflt);
  
  # % initial values for iteration
  # b=info.initb;
  
  deflt <- list(maxitr=20000,convitr=500,printitr=1000,conv=1e-6,initb=b0,maxaitr=-1)
  info  <- initstruct(info, deflt)
  
  b     <- info$initb
  
  # % initial estimate for a
  # maxaitr=info.maxaitr;
  maxaitr <- info$maxaitr
  
  # if(isfield(info,'inita'))
  #   a=info.inita;
  # else
  #   a=repmat(sqrt(sum(x.^2,2)/nt),1,nt);
  # if(maxaitr<0)
  #   maxaitr=10;
  # info.maxaitr=10;
  # end
  # end
  
  if (exists("inita", info)) {
    a              <- info$inita
  } else {
    a              <- pracma::repmat(sqrt(as.matrix(apply(x^2, 1, sum))/nt),1,nt)
    if (maxaitr < 0) {
      maxaitr      <- 10
      info$maxaitr <- 10
    }
  }

  # if(maxaitr>0)
  #   t=xw*b';
  #   for itr=1:maxaitr
  #       a=a.*t./(((a*b).*w2x)*b'+1e-100);
  # end
  # info.inita=a;
  # end

  if (maxaitr > 0) {
    t <- xw %*%  t(b)
    for (itr in 1:maxaitr) {
      a <- a * t / (((a %*% b) * w2x) %*% t(b) + 1e-100)
    }
    
    info$inita <- a
  }

  # % ensure iteration counts nest & convitr<=printitr<=maxitr
  # convitr=info.convitr;
  # printitr=info.printitr;
  # maxitr=info.maxitr;

  convitr  <- info$convitr
  printitr <- info$printitr
  maxitr   <- info$maxitr

  # if(printitr>maxitr)
  #   maxitr=ceil(maxitr/convitr)*convitr;
  # printitr=maxitr+1;
  # else
  #   printitr=ceil(printitr/convitr)*convitr;
  # maxitr=ceil(maxitr/printitr)*printitr;
  # end
  # info.printitr=printitr;
  # info.maxitr=maxitr;
  
  if (printitr > maxitr) {
    maxitr      <- pracma::ceil(maxitr/convitr)*convitr
    printitr    <- maxitr + 1
  } else {
    printitr    <- pracma::ceil(printitr/convitr)*convitr
    maxitr      <- pracma::ceil(maxitr/printitr)*printitr
  }
  
  info$printitr <- printitr
  info$maxitr   <- maxitr
  
  # % initial value rms
  # er=x-a*b;
  # err=[er./sdx;(b0-b)./sdb];
  # rms1=sqrt(sum(err(:).^2)/nerr);
  # conv=info.conv;
  
  er   <- x - a %*% b
  err  <- rbind(er / sdx, (b0-b) / sdb)
  rms1 <- sqrt(sum(err^2)/nerr)
  conv <- info$conv
  
  # % print heading & initial values
  # if(printitr<=maxitr)
  # rmsx=sqrt(sum(er(:).^2)/nx);
  # rmsxwt=sqrt(sum(sum((er./sdx).^2/nx)));
  # erb=b0-b;
  # rmsb=sqrt(sum(erb(:).^2)/nb0);
  # rmsbwt=sqrt(sum(sum((erb./sdb).^2))/nb0);
  # fprintf(['   itr    rmsx      rmsxwt     rmsb     '  ...
  #          'rmsbwt   drms       da         db\n'])
  # fprintf('%6i%#10.3g%#10.3g%#10.3g%#10.3g\n',  ...
  #         0,rmsx,rmsxwt,rmsb,rmsbwt)
  # end
  
  if (printitr<=maxitr) {
    rmsx   <- sqrt(sum(er^2)/nx)
    rmsxwt <- sqrt(sum((er/sdx)^2/nx))
    erb    <- b0 - b
    rmsb   <- sqrt(sum(erb^2)/nb0)
    rmsbwt <- sqrt(sum((erb/sdb)^2)/nb0)
    
    pracma::fprintf('   itr    rmsx      rmsxwt     rmsb     rmsbwt   drms       da         db\n')
    pracma::fprintf("%6i%#10.3g%#10.3g%#10.3g%#10.3g\n", 0,rmsx,rmsxwt,rmsb,rmsbwt)
    
  }
  
  # % main iteration loop
  # a1=a;
  # b1=b;
  
  a1 <- a
  b1 <- b
  
  for (itr in 1:maxitr) {
    # % update a & b
    # a=a.*(xw*b')./(((a*b).*w2x)*b'+1e-100);
    # b=b.*((a'*xw)+b0w)./(a'*((a*b).*w2x)+b.*w2b+1e-100);
    
    a <- a *(xw %*% t(b)) / (((a %*% b) * w2x)%*% t(b) + 1e-100)
    b <- b *((t(a) %*% xw) + b0w) / (t(a) %*% ((a %*% b) * w2x) + b * w2b+1e-100)
  
    # % check convergence occasionally
    # if(mod(itr,convitr)==0)
    #   daa=sqrt(sum((a1(:)-a(:)).^2)/na)/convitr;
    # dbb=sqrt(sum((b1(:)-b(:)).^2)/nb0)/convitr;
    # err=[(x-a*b)./sdx;(b0-b)./sdb];
    # rms=sqrt(sum(err(:).^2)/nerr);
    # drms=(rms1-rms)/convitr;
    
    if ((itr %% convitr) == 0) {
      daa  <- sqrt(sum((a1 - a)^2) / na) / convitr
      dbb  <- sqrt(sum((b1 - b)^2) / nb0) / convitr
      err  <- rbind((x - a %*% b) / sdx, (b0 - b) / sdb)
      rms  <- sqrt(sum(err^2) / nerr)
      drms <- (rms1 - rms ) / convitr
      
      # a=amatfactsd(x,sdx,b);
      # a1=a;
      # b1=b;
      # rms1=rms;
      # convind=abs(drms)<conv;
      
      # TODO: work on amatfactsd
      a       <- amatfactsd(x,sdx,b)
      a1      <- a
      b1      <- b
      rms1    <- rms
      convind <- abs(drms) < conv
      
      # % check print occasionally
      # if(mod(itr,printitr)==0 || convind || mod(itr,maxitr)==0)
      #   er=x-a*b;
      #   rmsx=sqrt(sum(er(:).^2)/nx);
      #   rmsxwt=sqrt(sum(sum((er./sdx).^2/nx)));
      #   erb=b0-b;
      #   rmsb=sqrt(sum(erb(:).^2)/nb0);
      #   rmsbwt=sqrt(sum(sum((erb./sdb).^2))/nb0);
      
      
      # if(printitr<=maxitr)
      #   fprintf(['%6i%#10.3g%#10.3g%#10.3g%#10.3g%',  ...
      #            '#11.3e%11.3e%#11.3e\n'],  ...
      #           itr,rmsx,rmsxwt,rmsb,rmsbwt,drms,daa,dbb)
      # end
      # end
      
      if ((itr %% printitr) == 0 || convind || (itr %% maxitr) == 0  ) {
        er     <- x - a %*% b
        rmsx   <- sqrt(sum(er^2) / nx)
        rmsxwt <- sqrt(sum((er / sdx)^2 / nx))
        erb    <- b0 - b
        rmsb   <- sqrt(sum(erb^2) / nb0)
        rmsbwt <- sqrt(sum((erb / sdb)^2) / nb0)
        
        if (printitr<=maxitr) {
          pracma::fprintf('   itr    rmsx      rmsxwt     rmsb     rmsbwt   drms       da         db\n')
          pracma::fprintf('%6i%#10.3g%#10.3g%#10.3g%#10.3g%#11.3e%11.3e%#11.3e\n', 
                          0,rmsx,rmsxwt,rmsb,rmsbwt,drms,daa,dbb)
        }
      }
      
      # % check convergence
      # if(convind)
      #   break
      
      # % check convergence
      if(convind) {
        break
      }
    }

  }
  
  # info.itr=itr;
  # info.conv=drms;
  # info.rms=rms;
  # info.rmsx=rmsx;
  # info.rmsxwt=rmsxwt;
  # info.rmsb=rmsb;
  # info.rmsbwt=rmsbwt;
  # 
  # return
  # end
  
  info$itr    <- itr
  info$conv   <- drms
  info$rms    <- rms
  info$rmsx   <- rmsx
  info$rmsxwt <- rmsxwt
  info$rmsb   <- rmsb
  info$rmsbwt <- rmsbwt
  
  
  result <- list(a=a,b=b,info=info)
  
  return(result)
  
  # ---- end ----  
}


