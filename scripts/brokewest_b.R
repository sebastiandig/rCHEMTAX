# ---- Everything up to here works as expected ----

# ---- set directory for saving  ----
root <- rprojroot::find_rstudio_root_file()
dir <- "/data/processed/"

# source scripts
source(paste0(root, "/scripts/chemtaxbrokewest.R"))
source(paste0(root,"/scripts/matfactuvw.R"))

# ---- extract data from chemtaxbrokewest ----
temp           <- chemtaxbrokewest(type = "norm") 
# df_matrix      <- temp$df_matrix      # pigment data
# df_matrix_sd   <- temp$df_matrix_sd   # df_matrix * 0.01 + 0.0003
# init_pig_ratio <- temp$init_pig_ratio # init_pig_ratio = ratio matrix
# pig_ratio_sd   <- temp$pig_ratio_sd   # init_pig_ratio * 0.1 w/ last col = 0.005
# taxa           <- temp$taxa           # name of taxa groups, comes from pigment ratio col 1
# pigm_sel       <- temp$pigm_sel       # pigment names, comes from pigment ratio colnames, keeps names where index is 1

s      <- temp$df      # pigment data
s_norm <- temp$df_norm   # df_matrix * 0.01 + 0.0003
f0     <- temp$pig_r_init # init_pig_ratio = ratio matrix
f_norm <- temp$pig_r_norm   # init_pig_ratio * 0.1 w/ last col = 0.005
taxa   <- temp$taxa           # name of taxa groups, comes from pigment ratio col 1
pigm   <- temp$pigm_sel 

# tic;[aaa,bbb,info]=matfactuvw(x2,b2,1,struct('maxitr',2000));toc
# %[aaa,bbb,info]=matfact5(x2,b2,struct('maxitr',200));


tictoc::tic()
temp <- matfactuvw(x = s_norm, b = f_norm, 1, info =list(maxitr = 2000))
aaa  <- temp$a
bbb  <- temp$b
info <- temp$info
tictoc::toc()


# t=x1(:,end)<0.3;
# x3=x2(t,:);

t  <- x1[,ncol(x1)] < 0.3
x3 <- x2[t,]

# semilogy(x1(:,end))
# title('CHEMTAXBROKEWest  Ch_a')
# xlabel('Sample number')
# ylabel('Chl-a')

x1.name <- colnames(x1)[ncol(x1)]

# allow y log style plotting from Matlab in semilogy
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

yticks = outer((1:10),(10^(-5:-1)))
xticks = outer((1:10),(10^(0:1)))


ggplot(x1, aes(y = !!sym(x1.name), x = seq(1, nrow(x1)))) +
  geom_point() +
  labs(title = 'CHEMTAXBROKEWest  Ch_a',
       x = 'Sample number',
       y ='Chl-a') +
  scale_y_log10(limits = c(1e-5, 1),
                labels = fancy_scientific,
                minor_breaks = yticks) +
  theme_bw()
