library("tictoc")
library("ggplot2")
library("pracma")

root <- rprojroot::find_rstudio_root_file()
source(paste0(root,"/scripts/permcalc.R"))
# TODO: matfactuvw
# source(paste0(root,"/scripts/matfactuvw.R"))

# bs=[1	0	0	1	1	0	1	1	1	1	1	0	1	1	0	0	1	0	1];
# bb={'chlc3'	'MgDVP'	'chlc2'	'chlc1'	'per'	'but'	'fuc'	'neox'	'prx'	'violax'	'hex'	'Mmal'	'alx'	'lut'	'dhlut'	'GyroxTotal'	'chl_b'	'np_chl_c2'	'chl_a'};
# b=[0	0.082815484	0	0	0	0	0	0.076017708	0.093775606	0.048498574	0	0.034069268	0	0.006383245	0.023810763	0	0.663053868	0	1
#    0	0.001674811	0	0	0	0	0	0.074384251	0	0.036350279	0	0	0	0.220905086	0	0	0.167089565	0	1
#    0	0.001280276	0.14910557	0	0	0	0	0	0	0	0	0	0.224584718	0	0	0	0	0	1
#    0	0.001108976	0.076997483	0.15	0	0	0.8	0	0	0	0	0	0	0	0	0	0	0	1
#    0.033	0	0.131	0	0	0	0.61	0	0	0	0	0	0	0	0	0	0	0	1
#    0	0.000974879	0.367132108	0	0.876791629	0	0	0	0	0	0	0	0	0	0	0	0	0	1
#    0.13	0.001284182	0.023	0	0	0.01	0.08	0	0	0	0.4	0	0	0	0	0	0	0.03	1
#    0.27	0.001120303	0.16	0	0	0.12	0.01	0	0	0	1.1	0	0	0	0	0	0	0.06	1];

# turns columns on/off (1/0) to be used in analysis
bs <- c(1, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1)
# bb <- c("chlc3", "MgDVP", "chlc2", "chlc1", "per", "but", "fuc", "neox", "prx", 
        # "violax", "hex", "Mmal", "alx", "lut", "dhlut", "GyroxTotal", "chl_b",
        # "np_chl_c2", "chl_a")

# reads pigment ratio start file and extracts col names
b <- readr::read_csv(paste0(root,"/scripts/b_pigment_ratios.csv"))
bb <- names(b)


# CHEMTAXBROKEWestx
# should just be a file to load
# i.e 
x <- readr::read_csv(paste0(root,"/scripts/CHEMTAXBROKEWestx.csv"))
xx <- colnames(x)

# xi=permcalc(bb(bs==1),xx);
# bi=find(bs);
# bb(bi)
# xx(xi)

# x1=x(:,xi);
# b1=b(:,bi);
# x2=x1./repmat(sum(x1,2),1,size(x1,2));
# b2=b1./repmat(sum(b1,2),1,size(b1,2));

xi <- permcalc(bb[bs == 1], xx)
bi <- which(bs == 1)
bb[bi]
xx[xi]

x1 <- x[,xi]
b1 <- b[,bi]
x2 <- x1/pracma::repmat(apply(x1, 1, sum),1,ncol(x1))
b2 <- b1/pracma::repmat(apply(b1, 1, sum),1,ncol(b1))



# tic;[aaa,bbb,info]=matfactuvw(x2,b2,1,struct('maxitr',2000));toc
# %[aaa,bbb,info]=matfact5(x2,b2,struct('maxitr',200));

tictoc::tic()
temp <- matfactuvw(x2, b2, 1, info = data.frame(maxitr = 2000))
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
 
ggplot(x1, aes(y = !!sym(x1.name), x = seq(1, nrow(x1)))) +
  geom_point() +
  labs(title = 'CHEMTAXBROKEWest  Ch_a',
       x = 'Sample number',
       y ='Chl-a') +
  theme_bw()

