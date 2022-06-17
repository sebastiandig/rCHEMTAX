permcalc <- function(a,b) {
################################################################################
#                                                                              # 
#               Calculate permutation to term a into b                         #
#                                                                              #    
################################################################################
# DESCRIPTION:
# --------
# This uses a vector of pigment names, `a`, that has been pre-filtered using an 
# index of 0/1 to select only the pigments wanted and uses these names to  
# identify the pigment names from the data set, `b`, that match. 
#
# If a pigment in `a` does not have a match in `b`, then a print message of  
# "{x name} is not found" will be displayed.
#
# Once the search is completed, the matching names are put into an output vector 
# `p` to be used on the data set for indexing.
#
# --------
# INPUTS:
# --------
# a       = vector of pigment names used that are filtered by an index variable 
#           of 0/1 (0 = not used, 1 = used) in the analysis
# b       = vector of column names from the data set used
# 
# --------
# OUTPUTS:
# --------
# p       = the selected names in `b` that matches `a`
# 
# --------
# NOTES:
# --------
# Original: 2010-02-09  Matlab7  W.Whiten
# --------
# References:
# --------
#
# --------
# Author:
# --------
# Sebastian Di Geronimo (Fri Jun 17 15:10:26 2022)
    
  
  # function p=permcalc(a,b)
  # % permcalc Calculate permutation to term a into b
  # %  2010-02-09  Matlab7  W.Whiten
  # 
  # n=length(a);
  # p=zeros(1,n);
  
  n <- length(a)
  p <- matrix(0, 1, n)

  # for i=1:length(a)
  # t=strmatch(a{i},b,'exact');
  # if(isempty(t))
  #   disp([a{i} '  not found'])
  # else
  #   p(i)=t;
  # end
  # end
  # 
  # return
  # end
  
  for (i in 1:n) {
    
    # t is the index in b that matches a
    t <- which(b == a[i])
    
    # if does not match, will message "is not found"
    if (identical(t,integer(0))) {
      
      message(paste(a[i], "is not found"))
      
    } else {
      
      # if matches will add index in b to p
      p[i] <- t
      
    }
  }
  
  return(p)
}