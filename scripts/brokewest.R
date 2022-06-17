# will need to update
# TODO: create this as a standalone function 


# % brokewest  Command file to fit sazchemtax090210 data
# %  2010-04-11  Matlab7  W.Whiten
# %
# % Input from chemtaxbrokewest function
# % Output to file brokewest.csv in comma separated format
# 
# % get data values
# [s,ssd,f0,fsd,taxa,pigm]=chemtaxbrokewest;
# 
# % fit the matrix factors
# [c,f,info]=nnmatfactsd(s,ssd,f0,fsd);
# 
# % scale the factors and original data
# [ss,cc,ff,rms]=normprod(s,c,f);
# 
# % write results to file brokewest.csv
# fid=fopen('brokewest.csv','w+');
# fprintf(fid,'chemtaxbrokewest\r\n\r\n');
# dataout(taxa,pigm,ff,fid);
# fprintf(fid,'\r\n');
# dataout(taxa,cc,fid);
# fclose(fid);
# 
# % % plot showing the effect of regularisation
# % regplot(s,ssd,f,fsd)
# %
# % % show converges from random starts for c
# % randstart(s,ssd,f0,fsd);
# % 
# % % bootstrap using parimetric log normal
# % bootln(s,ssd,f0,fsd);
# % 
# % % bootstrap non parametric on s
# % bootnp(s,ssd,f0,fsd)

# set directory for saving 
root <- rprojroot::find_rstudio_root_file()
dir <- "/data/processed/"

# TODO: source data
source(paste0(root, "scripts/chemtaxbrokewest.R"))
# source(paste0(root, "scripts/nnmatfactsd.R"))
# source(paste0(root, "scripts/normprod.R"))

# get data values
temp  <- chemtaxbrokewest() # should this be a function?
s     <- temp$s
ssd   <- temp$ssd
f0    <- temp$f0
fsd   <- temp$fsd
taxa  <- temp$taxa
pigm  <- temp$pigm

# fit the matrix factors
temp2 <- nnmatfactsd(s,ssd,f0,fsd)
c     <- temp2$c
f     <- temp2$f
info  <- temp2$info

# scale the factors and original data
temp3 <- normprod(s,c,f)
ss    <- temp$ss
cc    <- temp$cc
ff    <- temp$ff
rms   <- temp$rms

# # write results to file brokewest.csv
# fid=fopen('brokewest.csv','w+');
# fprintf(fid,'chemtaxbrokewest\r\n\r\n');
# dataout(taxa,pigm,ff,fid);
# fprintf(fid,'\r\n');
# dataout(taxa,cc,fid);
# fclose(fid)

# write results to file brokewest.csv
# TODO: check results of functions to be input to df1 and df2 
df1 <-  ff
colnames(df1) <- pigm
rownames(df1) <- tax

df2 <- cc
colnames(df2) <- tax

# Start a sink file with a CSV extension
sink(paste0(root, dir, 'brokewest.csv'))

# Write the first dataframe, with a title and final line separator
cat('chemtaxbrokewest\n\n')
write.csv(df1)
cat('\n')

# Write the 2nd dataframe to the same sink
write.csv(df2)

# Close the sink
sink()


do <- "n"
if (do == "y") {
  # TODO: sources
  # source(paste0(root, "scripts/regplot.R"))
  # source(paste0(root, "scripts/randstart.R"))
  source(paste0(root, "scripts/bootln.R"))
  source(paste0(root, "scripts/bootnp.R"))
  
  # plot showing the effect of regularisation
  regplot(s,ssd,f,fsd)
  
  # show converges from random starts for c
  randstart(s,ssd,f0,fsd);
   
  # bootstrap using parimetric log normal
  bootln(s,ssd,f0,fsd);
   
  # bootstrap non parametric on s
  bootnp(s,ssd,f0,fsd)
}

# examples of how to create similar brokewest
# pigm <- c("chlc3", "chlc1", "per", "fuc", "neox", "prx", "violax", "hex", "alx", "lut", "chl_b", "chl_a")
# df1 <- as.matrix(cbind(
#   c(0, 0, 0, 0, 0.07039, 0, 0.5443, 0.07962),
#   c(0, 0, 0, 0.1037, 0, 0, 0, 0),
#   c(0, 0, 0, 0, 0, 3.627, 0, 0),
#   c(0, 0, 0, 2.695, 1.964, 0, 0.1089, 0.0101),
#   c(0.1848, 4.94065645841247e-324, 0, 0, 0, 0, 0, 0),
#   c(0.1421, 0, 0, 0, 0, 0, 0, 0),
#   c(0.07833, 0.001081, 0, 0, 0, 0, 0, 0),
#   c(0, 0, 0, 0, 0, 0, 0.2006, 1.745),
#   c(0, 0, 0.4142, 0, 0, 0, 0, 0),
#   c(0.01048, 0.001551, 0, 0, 0, 0, 0, 0),
#   c(0.5926, 0.001043, 0, 0, 0, 0, 0, 0),
#   c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L)
# ))
# 
# 
# tax <- c("Prasinophytes", "Chlorophytes", "Cryptophytes", "Diatoms-A", "Diatoms-B", "Dinoflagellates-A", "Haptophytes-HiFe", "Haptophytes-LoFe")
# 
# colnames(df1) <- pigm
# rownames(df1) <- tax
# 
# 
# # Sample dataframes:
# pigm <- c("chlc3", "chlc1", "per", "fuc", "neox", "prx", "violax", "hex", "alx", "lut", "chl_b", "chl_a")
# df2 <- as.matrix(cbind(
#   V1 = c(0.002601,0.001525,0.002354,0.001152,
#          0.004978,0.03478,0.03137,0.04727,0.02407),
#   V2 = c(0.4948,0.5428,0.5386,0.492,0.4089,
#          0.297,0.2124,0.2302,0.2948),
#   V3 = c(0.009011,0.01054,0.004325,0.02006,
#          0.02126,0.06749,0.09231,0.07517,0.06479),
#   V4 = c(0.01309,0.01352,0.006997,0.02478,
#          0.02479,0.003574,0.003676,8.05e-05,0),
#   V5 = c(0.2894,0.2646,0.2791,0.2805,0.2532,
#          0.2341,0.2665,0.2181,0.1497),
#   V6 = c(0.01135,0.009738,0.01191,0.0104,
#          0.006588,0.001318,0.002852,0.001766,0.001105),
#   V7 = c(0.07932,0.06079,0.06654,0.06239,0.115,
#          0.1846,0.2606,0.2389,0.222),
#   V8 = c(0.09963,0.09632,0.0898,0.1088,0.1653,
#          0.1775,0.1308,0.1886,0.2435))
# )
# 
# colnames(df2) <- tax
# 
# # Start a sink file with a CSV extension
# sink('brokewest1.csv')
# 
# # Write the first dataframe, with a title and final line separator
# cat('chemtaxbrokewest\n\n')
# write.csv(df1)
# cat('\n')
# 
# 
# # Write the 2nd dataframe to the same sink
# write.csv(df2)
# 
# 
# # Close the sink
# sink()





