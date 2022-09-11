
library("ggplot2")

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
temp <- matfactuvw(.df = s_norm, .pig_r = f_norm, wb = 1, .info = list(maxitr = 20000))
aaa  <- temp$a
bbb  <- temp$b
info <- temp$info
tictoc::toc()


# t=x1(:,end)<0.3;
# x3=x2(t,:);

# ---- Everything up to here works as expected ----
s_last_col  <- s[,ncol(s)] < 0.3
x3 <- s_norm[s_last_col,]

# semilogy(x1(:,end))
# title('CHEMTAXBROKEWest  Ch_a')
# xlabel('Sample number')
# ylabel('Chl-a')

  # ---- initialize plot info ----
  source(paste0(root,"/scripts/fancy_scientific.R"))
  yticks = outer(1:10, 10^(-5:-1))
  xticks_minor = outer(1:10, 10^(0:1))

  
  s <- as.data.frame(s + 1)
  names(s) <-  pigm
  
for (i in seq(pigm)) {
  plt <- ggplot(s, aes_string(y = pigm[i], x = "seq(1, nrow(s))")) +
    
    geom_point() +
    labs(title = paste0('CHEMTAXBROKEWest (pigment = ', pigm[i],")"),
         x = 'Sample number',
         y =  bquote(ln(.(noquote(pigm[i]))~+~1))
         ) +
    scale_y_log10(limits = c(min(s[,i]), max(s[,i])),
                  labels = fancy_scientific,
                  minor_breaks = yticks) +
    theme_bw()
  print(plt)
}
  
