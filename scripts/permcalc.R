permcalc <- function(ratio_pig_filt, df_pig) {
################################################################################
#                                                                              # 
#                    Calculate permutation to turn a into b                    #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# This uses a vector of pigment names, `ratio_pig_filt`, that has been 
# pre-filtered using an index of 0/1 to select only the pigments wanted and uses 
# these names to  identify the pigment names from the data set, `df_pig`, that 
# match.
#
# If a pigment in `ratio_pig_filt` does not have a match in `df_pig`, then a 
# print message of "{x name} is not found" will be displayed. Also, a `stop` is
# initiated with a matrix of the selected pigments and which sample pigments 
# were not found.
#
# Once the search is completed, the matching names are put into an output vector 
# `new_idx` to be used on the data set for indexing.
#
# ---- INPUTS: -----------
# ratio_pig_filt  = vector of pigment names used that are filtered by an index  
#                   variable of 0/1 (0 = not used, 1 = used) in the analysis
# df_pig          = vector of column names from the data set used
# 
# ---- OUTPUTS: ----------
# new_idx         = the selected names in `df_pig` that matches `ratio_pig_filt`
# 
# ---- NOTES: ------------
# Original: 2010-02-09  Matlab7  W.Whiten
#
# ---- REFERENCES(s): ----
#
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (Fri Jun 17 15:10:26 2022)

  # ---- search for selected pigments in sample ----
  n       <- length(ratio_pig_filt) # number of pigments used in analysis
  new_idx <- matrix(0, 1, n) 
  
  for (i in seq(n)) {
    
    # t is the index in b that matches a
    match_idx    <- which(df_pig == ratio_pig_filt[i])
    
    # if does not match, will message "is not found"
    if (identical(match_idx,integer(0))) {

      message(paste(ratio_pig_filt[i], "is not found"))
      new_idx[i] <- NA
      
    } else {
      
      # if matches adds index in df_pig to new_idx
      new_idx[i] <- match_idx
      
    }
  }

  # ---- stop if any of the selected pigments were not found in the sample  ----
  if (!identical(ratio_pig_filt, df_pig[new_idx])) {
    message('Names do not match')
    print(cbind("Selected Pigments"= ratio_pig_filt,
                "Sample Pigments"  = df_pig[new_idx]))
    stop("Look at list of non-matching names!")
  }
  
  return(new_idx)
}