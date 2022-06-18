initstruct <- function(values,default) {
################################################################################
#                                                                              # 
#           Initialize struct with default values when not given              #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# Initialize the structure for matrix factorization 
# 
# --------
# INPUTS:
# --------
# values  = Struct with given values, should be a list or tibble?
# default = Struct with default values, should be a list or tibble?
#
# --------
# OUTPUTS:
# --------
# str     = Struct with default values added, should be a list or tibble?
#
# --------
# NOTES:
# --------
# Original: 2010-01-23  Matlab7  W.Whiten
#
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Sat Jun 18 17:26:29 2022) 
  

  # function str=initstruct(values,default)
  # % initstruct  Initiallise struct with default values when not given
  # %  2010-01-23  Matlab7  W.Whiten
  # %
  # % str=initstruct(value,default)
  # %  values  Struct with given values
  # %  default Struct with default values
  # %
  # %  str     Struct with default values added.
  # 
  
  # % copy default values
  # str=default;
  str <- default 
  
  # % copy new values to output
  # names=fieldnames(values);
  if (is.list(values)) {
    name <- names(values)
  } else {
    # if values is a tibble
    name <- values[,1]
  }
  
  # this adds the non-default info
  # for i=1:length(names)
  # str.(names{i})=values.(names{i});
  # end
  
  # this adds the non-default info
  for (i in 1:length(name)) {
    print(i)
    # str${{name}}
    str <- append(str, values[name[i]])
  }
  str
  # return
  # end
}

# testing parts
# values <- list('maxitr22'=111000)
# 
# deflt=list('maxitr'=1000,'printitr'=100,'conve'=1e-6,'conva'=1e-6,  
#              'convb'=1e-6)
# 
# tes <- initstruct(values, deflt)
# tes
