matfactuvw <- function(x,u=NULL,v=NULL,b,w=NULL,wb=1,info=NULL) {
################################################################################
#                                                                              # 
#               Non neg matrix factor with factored sd values                  #
#              Fast method in Matlab with approx standard deviations           #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# 
# --------
# INPUTS:
# --------
# x       = Matrix to be factored as a*b
# u       = Row sd factor for x (sd is u'*v)
# v       = Column sd factor for x
# b       = Stabilising value for b, non zero locations, & initial value
# w       = Row sd factor for for (sd is w'*v)
# wb      = Scalar weight factor for b errors (short call sequence only)
# info    = Information for calculation (to override default values)
#           .maxitr  Maximum number of iterations
#           .printitr Count of iterations between printing
#           .conva  Converge value for a
#           .convb  Converge value for b
#           .conve  Converge value for rms change
#           .inita  Initial a value
#           .initb  Initial b value
#
# --------
# OUTPUTS:
# --------
# a       = Left factor of x
# b       = Right factor of x
# info    = Initial info and
#           .itr  Iterations used
#           .conva Converge level for a
#           .convb Converge level for b
#           .conve Converge level for error
#
# --------
# NOTES:
# --------
# Original: 2010-02-21  Matlab7  W.Whiten
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Mon Jun 20 21:26:19 2022)

  
  # % matfactuvw  Non neg matrix factor with factored sd values
  # %              Fast method in Matlab with approx standard deviations
  # %  2010-02-21  Matlab7  W.Whiten
  # %
  # % (a,b,info]=matfactuvw(x,b,wb,info)
  # % [a,b,info]=matfactuvw(x,u,v,b,w,info)
  # %  x    Matrix to be factored as a*b
  # %  u    Row sd factor for x (sd is u'*v)
  # %  v    Column sd factor for x
  # %  b    Stabilising value for b, non zero locations, & initial value
  # %  w    Row sd factor for for (sd is w'*v)
  # %  wb   Scalar weight factor for b errors (short call sequence only)
  # %  info Information for calculation (to override default values)
  # %        .maxitr  Maximum number of iterations
  # %        .printitr Count of iterations between printing
  # %        .conva  Converge value for a
  # %        .convb  Converge value for b
  # %        .conve  Converge value for rms change
  # %        .inita  Initial a value
  # %        .initb  Initial b value
  # %
  # %  a    Left factor of x
  # %  b    Right factor of x
  # %  info Initial info and
  # %        .itr  Iterations used
  # %        .conva Converge level for a
  # %        .convb Converge level for b
  # %        .conve Converge level for error

  # % calc sd vectors for short calling sequence
  # if(nargin<5) - missing info
  #   if(nargin<4)
  #     info=struct();
  # else
  #   info=b;
  # end
  # b=u;
  # wb=v;
  # u=sqrt(sum(x.^2,2)/size(x,2));
  # v=0.01*sqrt(sum(x.^2,1)/size(x,1))+0.01;
  # w=wb*sqrt(sum(b.^2,2)/size(b,2));
  # 
  # % default for info
  # elseif(nargin<6)
  # info=struct();
  # end
  
  source(paste0(root,"/scripts/initstruct.R"))
  
  if (is.null(u)){
    # u=sqrt(sum(x.^2,2)/size(x,2));
    u <- sqrt(apply(x^2, 1, sum, na.rm=T) / dim(x)[2])
    
  }
  
  if (is.null(v)){
    # v=0.01*sqrt(sum(x.^2,1)/size(x,1))+0.01;
    v <- 0.01 * sqrt(apply(x^2, 2, sum, na.rm=T)/dim(x)[1]) + 0.01
  }
  
  if (is.null(w)){
    w <- wb * sqrt(apply(x^2, 1, sum, na.rm=T)/dim(b)[2]);
  }
  
  if (is.null(info)){
    info <- list()
  }
  
  # # % matrix sizes
  # [ns,np]=size(x);
  # [nt,nx]=size(b);
  # if(np~=nx); error('matfactuvw: matrix size error'); end
  
  ns <- dim(x)[1]
  np <- dim(x)[2]
  nt <- dim(b)[1]
  nx <- dim(b)[2]
  if(np != nx) stop('matfactuvw: matrix size error')
  
  # % adjust default values as given
  # deflt=struct('maxitr',1000,'printitr',100,'conve',1e-6,'conva',1e-6,  ...
  #              'convb',1e-6,'inita',ones(ns,nt)/nt,'initb',b);
  # info=initstruct(info,deflt);

  deflt <- list(maxitr=1000,printitr=100,conve=1e-6,conva=1e-6,
                convb=1e-6,inita=matrix(1,ns,nt)/nt,initb=b)
  info  <- initstruct(info, deflt)
  
  # % initail values for interation
  # aa=info.inita./repmat(u(:),1,nt);
  # 
  # bb=info.initb;
  # bb=bb./repmat(v(:)',nt,1);
  # bb=bb./repmat(v(end).*bb(:,end),1,np);
  # 
  # xuv=x./(u(:)*v(:)');
  # xb=[xuv;b./(w(:)*v(:)')];
  # sq=xb-[aa;diag(1./w)]*bb;
  # rms=sqrt(sum(sq(:).^2)/numel(sq));

  aa <- info$inita / pracma::repmat(u,1,nt)
  
  bb <- info$initb
  bb <- bb / pracma::repmat(t(v),nt,1)
  bb <- bb / pracma::repmat(v[nrow(v), ncol(v)] * bb[,ncol(bb)],1,np)

  xuv <- x / (matrix(u) %*% t(matrix(v)))
  xb <- rbind(xuv,
              b / (matrix(w) %*% t(matrix(v))))
  sq <- xb - rbind(aa,
                  pracma::Diag(1 / w)) %*% bb
  rms <- sqrt( sum(sq^2) / length(sq))

  
  # if(info.printitr<1e6)
  # rmsxuv=sqrt(sum(sum((xuv-aa*bb).^2))/numel(xuv));
  # rmsbb=sqrt(sum(sum((b./(w(:)*v(:)')-diag(1./w)*bb).^2))/  ...
  #       sum(b(:)~=0));
  # rmsx=sqrt(sum(sum((x-  ...
  #     (aa.*repmat(u,1,nt))*(bb.*repmat(v,nt,1))).^2))/numel(x));
  #
  # fprintf(['   itr     rmsx        rmswt       rmsb   '  ...
  #     'drms       daa        dbb\n'])
  # fprintf('%6i%#11.3g%#11.3g%#11.3g\n',0,rmsx,rmsxuv,rmsbb)
  #   end
    
  if (info$printitr < 1e6) {
    
    rmsxuv <- sqrt(sum((xuv - aa %*% bb)^2) / length(xuv))
    rmsbb  <-
      sqrt(sum((
        b / (matrix(w) %*%  t(matrix(v))) - pracma::Diag(1 / w) %*% bb) ^ 2) /
        sum(matrix(b)!= 0))
    
    rmsx   <-  
      sqrt(sum((
      x - (aa * pracma::repmat(u, 1, nt)) * (bb * pracma::repmat(v, nt, 1))) ^ 2) 
      / length(x))
    
    pracma::fprintf('   itr    rmsx      rmsxwt     rmsb     drms       daa       dbb\n')
    pracma::fprintf('%6i%#11.3g%#11.3g%#11.3g\n',0,rmsx,rmsxuv,rmsbb)
  }
  
  # for itr=1:info.maxitr
  # aa1=aa.*(xuv*bb')./(aa*(bb*bb'));
  # %aa1=aa.*(1+((xuv-aa*bb)*bb')./(aa*(bb*bb')));
  
  for (itr in 1:info$maxitr) {
    aa1 <- aa * (xuv %*% t(bb)) / (aa %*% (bb %*% t(bb)))
  
    # aw=[aa1;diag(1./w)];
    # bb1=bb.*(aw'*xb)./(aw'*(aw*bb));
    # %bb1=bb.*(1+(aw'*(xb-aw*bb))./(aw'*(aw*bb)));
    # %bb1=bb1./repmat(v(end).*bb(:,end),1,np);
    
    aw  <- rbind(aa1, pracma::Diag(1 / w))
    bb1 <- bb * (t(aw) %*% xb) / (t(aw) %*% (aw %*% bb))
    
    # daa=sqrt(sum((aa1(:)-aa(:)).^2)/numel(aa));
    # dbb=sqrt(sum((bb1(:)-bb(:)).^2)/numel(bb));
    # sq=xb-aw*bb1;
    # rms1=sqrt(sum(sq(:).^2)/numel(sq));
    # drms=rms-rms1;
    
    daa  <- sqrt(mean((aa1-aa)^2))
    dbb  <- sqrt(mean((bb1-bb)^2))
    sq   <- xb - aw %*% bb1
    rms1 <- sqrt(mean(sq^2))
    drms <- rms - rms1 
    
    # aa=aa1;
    # bb=bb1;
    # rms=rms1;
    
    aa  <- aa1
    bb  <- bb1
    rms <- rms1
    
    
    # if(mod(itr,info.printitr)==0)
    #   rmsxuv=sqrt(sum(sum((xuv-aa*bb).^2))/numel(xuv));
    # rmsbb=sqrt(sum(sum((b./(w(:)*v(:)')-diag(1./w)*bb).^2))/  ...
    #         sum(b(:)~=0));
        # rmsx=sqrt(sum(sum((x-  ...
        #     (aa.*repmat(u,1,nt))*(bb.*repmat(v,nt,1))).^2))/numel(x));
        # fprintf('%6i%#11.3g%#11.3g%#11.3g%#11.3e%11.3e%#11.3e\n',  ...
        #   itr,rmsx,rmsxuv,rmsbb,drms,daa,dbb)
        # end
        
        if(daa<info.conva||dbb<info.convb||drms<info.conve)
          break
        end
        end
    
        
    if (itr %% printitr == 0) {
      rmsxuv <- sqrt(sum((xuv - aa %*% bb)^2)/length(xuv))
      rmsbb <-
        sqrt(sum((
          b / (matrix(w) %*%  t(matrix(v))) - pracma::Diag(1 / w) %*% bb) ^ 2) /
            sum(matrix(b)!= 0))
      rmsx <-  
        sqrt(sum((
          x - (aa * pracma::repmat(u, 1, nt)) * (bb * pracma::repmat(v, nt, 1))) ^ 2) 
          / length(x))
      
      pracma::fprintf('%6i%#11.3g%#11.3g%#11.3g%#11.3e%11.3e%#11.3e\n',
              itr,rmsx,rmsxuv,rmsbb,drms,daa,dbb)
    }
    
    if(daa<info$conva || dbb<info$convb||drms<info$conve)
      break
        
    }
  
  
  # a=aa.*repmat(u,1,nt);
  # b=bb.*repmat(v,nt,1);
  # t=b(:,end);
  # b=b./repmat(t,1,np);
  # a=a.*repmat(t',ns,1);
  # 
  # info.itr=itr;
  # info.conva=daa;
  # info.convb=dbb;
  # info.conve=drms;
  
  a <- aa * pracma::repmat(u,1,nt)
  b <-bb * pracma::repmat(v,nt,1)
  t <-b[,ncol(b)]
  b <-b / pracma::repmat(t,1,np)
  a <-a *  pracma::repmat(t(t),ns,1)
  
  info$itr   <-itr
  info$conva <-daa
  info$convb <-dbb
  info$conve <-drms
        
  results <- list(a=a,b=b,t=t,info)     
  
  return(results)    
  # ---- end ----
}