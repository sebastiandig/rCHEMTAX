fancy_scientific <- function(l) {
################################################################################
#                                                                              # 
#                 Make better looking log based lables on graph                #
#                                                                              #    
################################################################################
# ---- DESCRIPTION: ------
# When given a vector of numbers, will change the form to be in scientific 
# notation as a string. This can be used in either x or y axis of a plot.
# 
# Example:
#         x <- c(0.1,1,10,100)
#         0.1   1.0  10.0 100.0
#         fancy_scientific(x)
#         expression('1'%*%10^-01, '1'%*%10^+00, '1'%*%10^+01, '1'%*%10^+02)
#
# The resulting form will look like: 1x10^-1, 1x10^0, 1x10^1, 1x10^2
# 
# May be use in ggplot as number labels on x or y axis marks
# scale_x_log10(limits = c(10^(df_xmin), 100),
#               labels = fancy_scientific,
#               breaks = 10^(df_xmin:2),
#               minor_breaks = xticks_minor,
#               expand = expansion(add = c(0, 0.15))
#               ) +
# 
# ---- INPUTS: -----------
# l    = vector of numbers to form expressions in the form of num x 10^y 
#
# ---- OUTPUTS: ----------
# text = expressions in scientific notation as a string
#
# ---- NOTES: ------------
#
# ---- REFERENCES(s): ----
# Original: https://stackoverflow.com/questions/11610377/how-do-i-change-the-formatting-of-numbers-on-an-axis-with-ggplot
# 
# ---- AUTHOR(s): --------
# Sebastian Di Geronimo (2022-08-31 17:25:11)
  
  # ---- turn vector in to character string in scientific notation ----
  l <- format(l, scientific = TRUE)
  
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  
  # return this as an expression
  parse(text = l)
}