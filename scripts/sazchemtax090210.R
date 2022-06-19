sazchemtax090210 <- function() {
################################################################################
#                                                                              # 
#                      Read phytoplankton data for saz                         #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# 
# --------
# INPUTS:
# --------
# NA
#
# --------
# OUTPUTS:
# --------
# s       = Matrix of samples by pigment readings
# ssd     = Standard deviations for s
# f0      = Initial matrix of taxa by pigment estimate 
# fsd     = Standard deviations for f
# taxa    = Cell array of taxa names
# pigm    = Cell array of pigment names
#
# --------
# NOTES:
# --------
# Original: 2010-03-21  Matlab7  W.Whiten
#   
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sat Jun 18 19:24:17 2022)
  
# function [s,ssd,f0,fsd,taxa,pigm]=sazchemtax090210
# % sazchemtax090210  Read phytoplanton data for saz
# %  2010-03-21  Matlab7  W.Whiten
# %
# % [s,ssd,f0,fsd,taxa,pigm]=sazchemtax090210
# %  Input data from this file and sazchemtax090210
# %
# %  s    Matrix of samples by pigment readings
# %  ssd  Standard deviations for s
# %  f0   Initial matrix of taxa by pigment estiments
# %  fsd  Standard deviations for f
# %  taxa Cell array of taxa names
# %  pigm Cell array of pigment names
# %
# % Pigments are selected according to the indicator array
# %  and origanised so columns of pigments correspond
# % For other data change function name and paste data in as bolow
# %  and similarly for a file of the sample data
# 
# % put sample data file name here
# % get s matrix and scol column names

  # library("readr")
  
  root <- rprojroot::find_rstudio_root_file()
  raw  <- "/data/raw/"
  
  # SAZS_CHEMTAX090210s
  s    <- as.matrix(read.csv(paste0(root, raw, "SAZS_CHEMTAX090210s.csv")))
  scol <- colnames(s)
  colnames(s) <- NULL
  
  # % SAZ_CHEMTAX090210b  Load chemtax data
  # ind=[1	0	0	1	1	1	1	1	1	1	1	1	1	1	1	1];
  # fcol={'Chl c3'	'MgDVP'	'Chl c2'	'Chl c1'	'Peridinin'	'ButFuc'	'Fuc'	'pras'	'violax'	'Hex'	'allox'	'Zea'	'Lut'	'Chlb'	'np-chl c2'	'chl_a'};
  # 
  # f0=[0	0.062	0	0	0	0	0	0.245	0.054	0	0	0.058	0.021	0.704	0	1
  #     0	0.001674811	0	0	0	0	0	0	0.049	0	0	0.032	0.171	0.315	0	1
  #     0	0	0	0	0	0	0	0	0	0	0	0	0	0.377	0	1
  #     0	0.001280276	0.196	0	0	0	0	0	0	0	0.379	0	0	0	0	1
  #     0	0.001108976	0.045	0.023	0	0	0.623	0	0	0	0	0	0	0	0	1
  #     0.083	0	0.284	0	0	0	0.998	0	0	0	0	0	0	0	0	1
  #     0	0.000974879	0.218	0	0.558	0	0	0	0	0	0	0	0	0	0	1
  #     0.177	0.009	0.209	0	0	0.005	0.229	0	0	0.47	0	0	0	0	0.087	1
  #     0.171	0.029	0.192	0	0	0.103	0.3	0	0	0.371	0	0	0	0	0.058	1
  #     0.035	0.006	0.093	0	0	0.061	0.19	0	0	0.188	0	0	0	0	0	1
  #     0	0	0	0	0	0	0	0	0	0	0	0.636	0	0	0	1];
  
  # % b = [0	0.080440439	0	0	0	0	0	0.321291305	0.062825477	0	0	0.067087644	0.018587707	0.63519521	0	1
  #        % 0	0.001924497	0	0	0	0	0	0	0.03596378	0	0	0.025838928	0.224682924	0.278000009	0	1
  #        % 0	0	0	0	0	0	0	0	0	0	0	0	0	0.404873771	0	1
  #        % 0	0.001093119	0.143770289	0	0	0	0	0	0	0	0.400559428	0	0	0	0	1
  #        % 0	0.000734897	0.045799117	0.022933037	0	0	0.537369433	0	0	0	0	0	0	0	0	1
  #        % 0.059114731	0	0.311679121	0	0	0	1.107005575	0	0	0	0	0	0	0	0	1
  #        % 0	0.000976137	0.195690774	0	0.420229494	0	0	0	0	0	0	0	0	0	0	1
  #        % 0.225790245	0.009532848	0.215835286	0	0	0.006200425	0.18661012	0	0	0.335536215	0	0	0	0	0.069178843	1
  #        % 0.214704353	0.021671631	0.170104421	0	0	0.133045085	0.296671195	0	0	0.392466057	0	0	0	0	0.052025435	1
  #        % 0.028987712	0.005547827	0.109793135	0	0	0.051489242	0.210102747	0	0	0.224796597	0	0	0	0	0	1
  #        % 0	0	0	0	0	0	0	0	0	0	0	0.774158857	0	0	0	1];
  # 
  # taxa={'Prasinophytes'
  #   'Chlorophytes'
  #   'Euglenophytes'
  #   'Cryptophytes'
  #   'Diatoms-A'
  #   'Diatoms-B'
  #   'Dinoflag-A'
  #   'Hapto-A'
  #   'Hapto-B'
  #   'Dinoflag-B'
  #   'Cyanobact'};
  
  # read pigment ratios, extract phyto names, pigment names, and ratios as
  # separate variables
  ind <- c(1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
  temp <- read.csv(paste0(root,"/scripts/saz_pigment_ratios.csv"))
  f0 <- as.matrix(temp[,-1])
  fcol <- colnames(f0)
  taxa <- temp[,1]
  colnames(f0) <- NULL
  rm(temp)


  source(paste0(root,"/scripts/permcalc.R"))
  # % find columns that match
  # si=permcalc(fcol(ind==1),scol);
  # fi=find(ind);
  si <- permcalc(fcol[ind == 1], scol)
  fi <- which(ind == 1)
  
  # if(~isequal(fcol(fi),scol(si)))
  #   disp('Names do not match')
  # disp([fcol(fi),scol(si)])
  # return
  # end
  # pigm=fcol(fi);
  
  if (!all.equal(fcol[fi],scol[si])) {
    message('Names do not match')
    print(cbind(fcol[fi],scol[si]))
    # TODO: add exit?
    # return
  }
  
  pigm <- fcol[fi]
  
  
  # % get columns that match
  # s=s(:,si); %#ok<NODEF>
  # f0=f0(:,fi);
  
  s  <- s[, si]
  f0 <- f0[, fi]
  
  # % set sd values
  # ssd=s*0.01+0.0003;
  # fsd=f0*0.1;
  # fsd(:,end)=0.005;
 
  ssd <- s * 0.01 + 0.0003
  fsd <- f0*0.1
  fsd[,ncol(fsd)] <- 0.005
  
  # return
  # end
  
  result <-  list(s=s,ssd=ssd,f0=f0,fsd=fsd,taxa=taxa,pigm=pigm)
  
}