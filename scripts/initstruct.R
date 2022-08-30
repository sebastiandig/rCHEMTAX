initstruct <- function(values, default) {
################################################################################
#                                                                              # 
#             Initialize options with default values when not given            #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# Initialize the structure for matrix factorization.
# Default values:
#                maxitr   = 20000;      maximum iterations
#                convitr  = 500;        convergence iterations?
#                printitr = 1000,       print results every 1000 iterations
#                conv     = 1e-6,       convergence tolerance
#                initb    = pig_r_init, initial pigment ratios matrix
#                maxaitr  = -1          maximum for initial iterations for a
# 
# ---- INPUTS: -----------
# values  = Struct with given values, should be a list
# default = Struct with default values, should be a list
#
# ---- OUTPUTS: ----------
# str     = Struct with default values added, should be a list or tibble?
#
# ---- NOTES: ------------
# Original: 2010-01-23  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Sat Jun 18 17:26:29 2022) 

  # copy default values
  str <- default 
  
  # copy new values to output
  if (is.list(values)) {
    name <- names(values)
  } else {
    stop("values needs to be a list")
  }
  
  # this adds the non-default info
  for (i in 1:length(name)) {
    str[ {{name[i]}} ] <- values[name[i]]
  }
  
  # return structure
  str

}

