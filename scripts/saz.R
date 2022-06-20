# will need to update
# TODO: create this as a standalone function 

# % saz Command file to fit sazchemtax090210 data
# %  2010-04-11  Matlab7  W.Whiten
# %
# % Input from sazchemtax090210 function
# % Output to file saz.csv in comma separated format
# 
# % get data values
# [s,ssd,f0,fsd,taxa,pigm]=sazchemtax090210;
# 
# % fit the matrix factors
# [c,f,info]=nnmatfactsd(s,ssd,f0,fsd);
# 
# % scale the factors and original data
# [ss,cc,ff,rms]=normprod(s,c,f);
# 
# % write result to file saz.csv
# fid=fopen('saz.csv','w+');
# fprintf(fid,'sazchemtax090210\r\n\r\n');
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
source(paste0(root, "/scripts/sazchemtax090210.R"))
source(paste0(root, "/scripts/nnmatfactsd.R"))
source(paste0(root, "/scripts/normprod.R"))

# get data values
temp  <- sazchemtax090210() # should this be a function?
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

# write results to file brokewest.csv
# TODO: check results of functions to be input to df1 and df2 
df1 <-  ff
colnames(df1) <- pigm
rownames(df1) <- tax

df2 <- cc
colnames(df2) <- tax

# Start a sink file with a CSV extension
sink(paste0(root, dir, 'saz.csv'))

# Write the first dataframe, with a title and final line separator
cat('sazchemtax090210\n\n')
write.csv(df1)
cat('\n')

# Write the 2nd dataframe to the same sink
write.csv(df2)

# Close the sink
sink()


do <- "n"
if (do == "y") {
  # TODO: sources
  source(paste0(root, "/scripts/regplot.R"))
  source(paste0(root, "/scripts/randstart.R"))
  source(paste0(root, "/scripts/bootln.R"))
  source(paste0(root, "/scripts/bootnp.R"))
  
  # plot showing the effect of regularisation
  regplot(s,ssd,f,fsd)
  
  # show converges from random starts for c
  randstart(s,ssd,f0,fsd);
  
  # bootstrap using parimetric log normal
  bootln(s,ssd,f0,fsd);
  
  # bootstrap non parametric on s
  bootnp(s,ssd,f0,fsd)
}