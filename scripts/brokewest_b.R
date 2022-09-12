
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

# ---- run fast method of factor analysis  ----
tictoc::tic()
temp <- matfactuvw(.df = s_norm, .pig_r = f_norm, wb = 1, .info = list(maxitr = 2000))
aaa  <- temp$a
bbb  <- temp$b
info <- temp$info
chlor_norm <- temp$chlor_norm
tictoc::toc()

# ---- plotting? ----
s_last_col  <- s[,ncol(s)] < 0.3
x3 <- s_norm[s_last_col,]

# ---- initialize plot info ----
source(paste0(root,"/scripts/fancy_scientific.R"))
yticks = outer(1:10, 10^(-5:-1))
xticks_minor = outer(1:10, 10^(0:1))


# ---- plot each pigment as ln(pigment + 1) ----
s <- as.data.frame(s)
names(s) <-  pigm
s1 <- s + 1
  
for (i in seq(pigm)) {
  plt <- ggplot(s1, aes_string(y = pigm[i], x = "seq(1, nrow(s))")) +
    
    geom_point() +
    labs(title = paste0('CHEMTAXBROKEWest (pigment = ', pigm[i],")"),
         x = 'Sample number',
         y =  bquote(ln(.(noquote(pigm[i]))~+~1))
         ) +
    scale_y_log10(limits = c(min(s1[,i]), max(s1[,i])),
                  labels = fancy_scientific,
                  minor_breaks = yticks) +
    theme_bw()
  print(plt)
}
  
# ---- plot chlor-a less than 0.3 ug/L ----
s_last_col  <- s[,ncol(s)] < 0.3
s2          <- s[s_last_col,]

ggplot(s2, aes(x = seq(1, nrow(s2)), y = chl_a)) +
  geom_point() +
  labs(
    title = paste0('CHEMTAXBROKEWest (pigment = ', pigm[length(pigm)], ")"),
    x     = 'Sample number',
    y     =  bquote(ln(.(noquote(pigm[length(pigm)])) ~ + ~ 1))
) +
  scale_y_log10(
    limits       = c(min(s2$chl_a), max(s2$chl_a)),
    labels       = fancy_scientific,
    minor_breaks = yticks
  ) +
  theme_bw()

# ---- plot chlor and predicated chlor ----
# TODO: figure out how to convert to s
pig_est   <- aaa %*% bbb * s$chl_a
chlor_est <- pig_est[, ncol(pig_est)]

df_est    <- as.data.frame(cbind(s, chlor_est))

ggplot() +
  geom_point(data   = df_est, aes(x = chl_a, y = chlor_est)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")


# ---- plot noralized s and s ----
x3 <- s_norm[s_last_col,ncol(s_norm)] # not used for anything?
s1_norm <- as.data.frame(cbind(s2, x3))

ggplot() +
  geom_point(data = s1_norm, aes(x = chl_a, y = x3))

