chemtaxbrokewest <- function() {
################################################################################
#                                                                              # 
#               Read phytoplankton data for chemtaxbrokewest                    #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
#
# --------
# INPUTS:
# --------
# Nothing, all is set
# --------
# OUTPUTS:
# --------
# s       = Matrix of samples by pigment readings
# ssd     = Standard deviations for s
# f0      = Initial matrix of taxa by pigment estiments
# fsd     = Standard deviations for f
# taxa    = Cell array of taxa names
# pigm    = Cell array of pigment names
#
# --------
# NOTES:
# --------
# Original: 2010-03-21  Matlab7  W.Whiten
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Mon Jun 13 18:21:17 2022)
 
# function [s,ssd,f0,fsd,taxa,pigm]=chemtaxbrokewest
# % chemtxbrokewest  Read phytoplanton data for chemtaxbrokewest
# %  2010-03-21  Matlab7  W.Whiten
# %
# % [s,ssd,f0,fsd,taxa,pigm]=chemtaxbrokewest
# %  Input data from this filw and CHEMTAXBROKEWests
# %
# %  s    Matrix of samples by pigment readings
# %  ssd  Standard deviations for s
# %  f0    Initial matrix of taxa by pigment estiments
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
  
  root <- rprojroot::find_rstudio_root_file()
  raw <- "/data/processed/"
  
# CHEMTAXBROKEWests # Pigment concentration in samples
  # should make this a function?
  s <- read.csv(paste0(root, raw, "CHEMTAXBROKEWests.csv"))
  scol <- colnames(s)
  
# % ind tell which columns can be used
# ind=[1	0	0	1	1	0	1	1	1	1	1	0	1	1	0	0	1	0	1];
 ind <- c(1L, 0L, 0L, 1L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 1L, 0L, 0L, 1L, 0L, 1L)


# % column names of pigments for f matrix
# fcol={'chlc3'	'MgDVP'	'chlc2'	'chlc1'	'per'	'but'	'fuc'	'neox'	'prx'	'violax'	'hex'	'Mmal'	'alx'	'lut'	'dhlut'	'GyroxTotal'	'chl_b'	'np_chl_c2'	'chl_a'};
  fcol <- c("chlc3",
            "MgDVP",
            "chlc2",
            "chlc1",
            "per",
            "but",
            "fuc",
            "neox",
            "prx",
            "violax",
            "hex",
            "Mmal",
            "alx",
            "lut",
            "dhlut",
            "GyroxTotal",
            "chl_b",
            "np_chl_c2",
            "chl_a"
            )

# % f matrix taxa by pigment
# f0=[0	0.082815484	0	0	0	0	0	0.076017708	0.093775606	0.048498574	0	0.034069268	0	0.006383245	0.023810763	0	0.663053868	0	1
#     0	0.001674811	0	0	0	0	0	0.074384251	0	0.036350279	0	0	0	0.220905086	0	0	0.167089565	0	1
#     0	0.001280276	0.14910557	0	0	0	0	0	0	0	0	0	0.224584718	0	0	0	0	0	1
#     0	0.001108976	0.076997483	0.15	0	0	0.8	0	0	0	0	0	0	0	0	0	0	0	1
#     0.033	0	0.131	0	0	0	0.61	0	0	0	0	0	0	0	0	0	0	0	1
#     0	0.000974879	0.367132108	0	0.876791629	0	0	0	0	0	0	0	0	0	0	0	0	0	1
#     0.13	0.001284182	0.023	0	0	0.01	0.08	0	0	0	0.4	0	0	0	0	0	0	0.03	1
#     0.27	0.001120303	0.16	0	0	0.12	0.01	0	0	0	1.1	0	0	0	0	0	0	0.06	1];
  # I created a csv for this instead
  f0 <- readr::read_csv(paste0(root,"/scripts/pigment_ratios.csv"), col_names = F)


# % row names for f matrix
# taxa={'Prasinophytes'
#   'Chlorophytes'
#   'Cryptophytes'
#   'Diatoms-A'
#   'Diatoms-B'
#   'Dinoflagellates-A'
#   'Haptophytes-HiFe'
#   'Haptophytes-LoFe'};
  # TODO: this could be a csv as well
  taxa <- c('Prasinophytes',
    'Chlorophytes',
    'Cryptophytes',
    'Diatoms-A',
    'Diatoms-B',
    'Dinoflagellates-A',
    'Haptophytes-HiFe',
    'Haptophytes-LoFe')

  # source(paste0(root,"scripts/permcalc.R"))
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
  
# % get columns that match
# s=s(:,si); %#ok<NODEF>
# f0=f0(:,fi);
  s  <- s[, si]
  f0 <- f0[, fi]
  
# % set sd values
# ssd=s*0.01+0.0003;
# fsd=f0*0.1;
# fsd(:,end)=0.005;
  ssd <- s*0.01+0.0003
  fsd <- f0*0.1
  fsd[,ncol(fsd)] <- 0.005
  
# return
# end
# 
   
 result <-  list(s,ssd,f0,fsd,taxa,pigm)

}

