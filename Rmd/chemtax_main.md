CHEMTAX Main
================
Sebastian DiGeronimo
6/2/2022

### Intro:

This will serve as my starting point for converting CHEMTAX V2.0 into R.

\#TODO: fix this: [CHEMTAX
Document](./Re__CHEMTAX_Update_%5BSEC=OFFICIAL%5D/Re__CHEMTAX_Update_%5BSEC=OFFICIAL%5D/CHEMTAX%202%20Octave%20version.docx)

\# Files that need to be converted double tilde \~ will show done

## Formulas:

Iterations without standard deviations:

$$
a = \frac {a \times(x \cdot b^T)}
          {a \cdot(b \cdot b^T) + 10^{-100}} 
\ ; \ 
b = \frac{b \times (a^T \cdot x + b_0)}
         {(a^T \cdot(a \cdot b) + b) + 10^{-100}}
$$

Iterations with standard deviations:

$$
\begin{aligned}
& a_{new} = \frac {a \times(x_{w} \cdot b^T)}
                {((a \cdot b ) \times x_{\sigma^2}^{-1}) \cdot b^T + 10^{-100}} 
\ ; \ 
b_{new} = \frac {b \times (a^T \cdot x_w) + b_{0w}} 
                {a^T \cdot ((a \cdot b) \times x_{\sigma^2}^{-1}) + b \times b_{\sigma^2}^{-1} + 10^{-100}} \\\\
& \text {where}; \\
& x = \text{pigment concentations for each sample}; \\
& x_{\sigma^2}^{-1} = \frac {1}{x_{\sigma}^2} = \text{inverse-variance weighting of }\ x; \\ 
& x_w = x \times x_{\sigma^2}^{-1} = \text{weighted } x \text{; } \\
& a = \text{taxa contribution}; \\
& b = \text{pigment ratio matrix};\\ 
& b_{\sigma^2}^{-1} = \frac {1}{b_{0\sigma}^2} = \text{inverse-variance weighting of initial}\ b; \\ 
& b_{0w} = b_0 \times b_{\sigma^2}^{-1}  = \text{weighted initial}\ b \\
\end{aligned}
$$

Non-negative Least Squares:

- Minimize $a$ matrix

$$
\begin{aligned}
& min||x_i - \hat{a_i} \cdot b_i|| \\
& i = \text{sample row}; \\
& x_i = (\frac {x_i} {x_{\sigma_i}})^2; \\
& b_i = \frac {b^T} {x_{\sigma_i}}; \\
& \hat{a_i} = \text{estimated taxa contribution per sample}\\
\end{aligned}
$$

- Minimize $b$ matrix

  $$
  \begin{aligned}
  & min||x_{ij} - a_{ik} \cdot \hat{b_{jk}}|| \\
  & i = \text{sample};\
  j = \text{pigment};\
  k = \text{taxa}; \\
  & k,j \in \Phi; \ 
  \Phi = \text{set of}\ k,j\ \text{non zero elements in pigment ratio matrix} \\
  & x_{j} = 
  \begin{pmatrix} 
          \frac {x_{ij}} {x_{\sigma_{ij}}} \\
          \frac {b_{jk}} {b_{\sigma_{jk}}}
  \end{pmatrix} ;\\
  & a_{kj} = 
  \begin{pmatrix}  
          \frac {a_{ijk}} {x_{\sigma_{ij}}} \\
          diag( b_{\sigma_{kj}}^{-1} )
  \end{pmatrix} ;\\
  & \hat {b_{jk}} = \text {estimated pigment ratio matrix} \\
  \end{aligned}
  $$

## List of functions or files

List of ones to go through again: bootln bootnp CHEMTAXBROKEWestb
matfactuvw randstart regplot saz sazchemtax090210

21 - not all are .m files with functions o 0 mostly done \* 20 finished
! 0 not started x 2 not needed

|                              Function                               | Status                                      | Connections                                                                                                         | Notes                                                                                       |
|:-------------------------------------------------------------------:|---------------------------------------------|---------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|
|                              *adddir*                               | <span style="color:red">**NA**</span>       |                                                                                                                     | I don’t think I need this                                                                   |
|                            *amatfactsd*                             | <span style="color:green">**Done**</span>   |                                                                                                                     | solve for a                                                                                 |
|                            *bmatfactsd*                             | <span style="color:green">**Done**</span>   |                                                                                                                     | solves for b                                                                                |
|                              *bootln*                               | <span style="color:yellow">**Works**</span> | *nnmatfactsd*                                                                                                       | close to finishing                                                                          |
|                              *bootnp*                               | <span style="color:yellow">**Works**</span> | *nnmatfactsd*                                                                                                       | close to finishing                                                                          |
|                             *brokewest*                             | <span style="color:green">**Done**</span>   | *chemtaxbrokewest*,                                                                                                 |                                                                                             |
| *nnmatfactsd, normprod*, *regplot*, *randstart*, *bootln*, *bootnp* |                                             |                                                                                                                     |                                                                                             |
|                         *chemtaxbrokewest*                          | <span style="color:green">**Done**</span>   | *permcalc,* file for *s* and *scol* as *CHEMTAXROKEWests*, file for *f0*, *taxa* and *fcol* as *pigment_ratios.csv* | calls data and pigment ratio matrix, selects pigments for both and creates std dev matrices |
|                         *CHEMTAXBROKEWestb*                         | <span style="color:green">**Done**</span>   | *matfactuvw*, *permcalc*, fo;e fpr *b* as brokewest_b\_*pigment_ratio.csv*, file for *x* as *CHEMTAXBROKEWestx.csv* | TODO: understand output                                                                     |
|                         *CHEMTAXBROKEWests*                         | <span style="color:green">**Done**</span>   |                                                                                                                     | as *.csv* in *data/raw/*                                                                    |
|                         *CHEMTAXBROKEWestx*                         | <span style="color:green">**Done**</span>   |                                                                                                                     | as *.csv* in *data/raw/*                                                                    |
|                              *dataout*                              | <span style="color:red">**NA**</span>       |                                                                                                                     | probably don’t need                                                                         |
|                              Formulas                               | <span style="color:green">**Done**</span>   |                                                                                                                     | doesn’t connect to anything and is not a function, not sure what do do with it              |
|                            *initstruct*                             | <span style="color:green">**Done**</span>   |                                                                                                                     | creates options variable for factorization                                                  |
|                            *matfactuvw*                             | <span style="color:green">**Done**</span>   | *initstruct*                                                                                                        |                                                                                             |
|                            *nnmatfactsd*                            | <span style="color:green">**Done**</span>   | *amatfactsd*, *initstruct*                                                                                          | <span style="color:green">**MAIN script for factor analysis**</span>                        |
|                             *normprod*                              | <span style="color:green">**Done**</span>   |                                                                                                                     | normalizes to selected column                                                               |
|                             *permcalc*                              | <span style="color:green">**Done**</span>   |                                                                                                                     |                                                                                             |
|                             *randstart*                             | <span style="color:green">**Done**</span>   | *nnmatfactsd*                                                                                                       |                                                                                             |
|                              *regplot*                              | <span style="color:green">**Done**</span>   |                                                                                                                     |                                                                                             |
|                                 saz                                 |                                             | *sazchemtax090210* *nnmatfactsd* *normprod* *regplot* *randstart* *bootln* *bootnp*                                 |                                                                                             |
|                         *sazchemtax090210*                          |                                             | as *.csv* in *data/raw/*                                                                                            |                                                                                             |

## To test `nnmatfactsd` and `amatfactsd` - start with

    - brokewest.R to line 57, that'll be the inputs to `nnmatfactsd`

# Matlab to R notes:

- `length(x)` = largest of the dimensions, so if 100 x 4 = 100, or 4 x
  100 = 100
- `*` is matrix multiply,
- `.*` & `./` are matrix element by element multiply & divide operations
- can use `pracma::rand()` for rand
- `repmat` is used to replicate a matrix by `x * y`, ie
  `repmat({matrix 5 x 4}, 2, 4)` result would be matrix 2 x 16, so 8
  replicates of the original matrix where stacked twice and multiplied
  by 4
- `mod(a,m)` is same as `a %% m` in R
- `df(:)` puts matrix into 1 column
- `df(end)` gets the last digit, bottom right of matrix

## Normalization to Chlor-a

I think this is done in `matfactuvw` using `pracma::Diag(1/chlor-a col)`
then the dot product with df

\#TODO: something to include is: This would load internal ratios file so
that you have a basic pigment ratio var
`pigment_ratio <- system.file("", package = "rCHEMTAX")`

## Plot in log notation on axes

\# allow y or x log style plotting from Matlab in semilogy, semilogx,
loglog

$$
\text{i.e. }x\times10^{n}
$$

``` r
fancy_scientific <- function(l) { 
# turn in to character string in scientific notation 

# quote the part before the exponent to keep all the digits 
l <- format(l, scientific = TRUE) 

# turn the 'e+' into plotmath format 
l <- gsub("(.)e", "'\1'e", l) 
l <- gsub("e", "%%10", l) 

# return this as an expression 
parse(text=l) }

# also, needs to set tick marks for major/minor; here 10e-5 to 10e10
yticks = outer((1:10),(10^(-5:-1))) 
xticks = outer((1:10),(10^(0:1)))

# include:
semi-logx scale_x_log10(limits = c(1, 100), labels = fancy_scientific, minor_breaks = xticks) + 
semi-logy scale_y_log10(limits = c(1e-5, 1), labels = fancy_scientific, minor_breaks = yticks) +
```

# Functions from `pracma`

- `repmat`
- `fprintf`
- `ceil` - not anymore, was alias for `base::ceiling`
- `rand`
- `Diag` - not anymore
- `lsqnonneg` - could replace with `limsolve`

## Testing stuff/setup project

``` r
# will create a set directory if does not exists
# useful for new projects
mainDir <- rprojroot::find_rstudio_root_file()
subDir <-
    c("data/raw",
      "data/processed",
      "data/plots",
      "data/metadata",
      "Rmd",
      "scripts")

fs::dir_create(path = paste0(mainDir,"/",subDir))
rm(mainDir, subDir)
```

# Normalize by last column

when normalizing by a column, if you take last column and create into
square matrix with diagonal being each value in last column and take the
inverse, you can take the dot product of the diagnoal matrix and the
original dataframe to get a normalized matrix by last column

Ex: dataset = 6x12 Take last column to normalize to diag(1 / col 12)
%\*% dataset will results if took last column and divided by each value
in the same row

``` r
test_norm <- pracma::rand(6, 12)

# last column should be a 1 because normalized to last column
# 2 methods:
# diagonol approach
diag(1/test_norm[,ncol(test_norm)]) %*% test_norm
```

    ##           [,1]      [,2]     [,3]      [,4]     [,5]      [,6]      [,7]
    ## [1,] 1.3509715 0.6427029 0.099324 0.4372269 0.474389 1.1037331 0.5238595
    ## [2,] 2.6069082 0.6543306 3.560185 3.1979000 1.524747 0.5186843 3.1476803
    ## [3,] 7.1320460 3.0969873 1.339884 5.9831552 5.742415 3.5512205 5.0871473
    ## [4,] 0.5572088 1.6275385 1.698656 0.1816916 1.550892 1.2916801 1.4549099
    ## [5,] 2.0082564 2.9530885 3.551854 1.0089437 3.294445 1.2265588 2.4998624
    ## [6,] 1.0478676 0.1059490 1.116854 0.8101214 0.694529 1.0466708 0.3551375
    ##          [,8]      [,9]     [,10]     [,11] [,12]
    ## [1,] 1.087794 0.8277188 1.0779756 0.5651250     1
    ## [2,] 1.072401 0.7409449 0.4494975 2.8441674     1
    ## [3,] 4.888614 0.2226610 6.6084202 4.6273951     1
    ## [4,] 1.060282 0.1368132 0.1094236 1.1665615     1
    ## [5,] 1.644214 2.2686004 2.8051097 1.9136866     1
    ## [6,] 1.056410 0.3154483 0.7324739 0.5862605     1

``` r
# apply method *more confusing
t(apply(test_norm, 1, function(x) x / x[length(x)]))
```

    ##           [,1]      [,2]     [,3]      [,4]     [,5]      [,6]      [,7]
    ## [1,] 1.3509715 0.6427029 0.099324 0.4372269 0.474389 1.1037331 0.5238595
    ## [2,] 2.6069082 0.6543306 3.560185 3.1979000 1.524747 0.5186843 3.1476803
    ## [3,] 7.1320460 3.0969873 1.339884 5.9831552 5.742415 3.5512205 5.0871473
    ## [4,] 0.5572088 1.6275385 1.698656 0.1816916 1.550892 1.2916801 1.4549099
    ## [5,] 2.0082564 2.9530885 3.551854 1.0089437 3.294445 1.2265588 2.4998624
    ## [6,] 1.0478676 0.1059490 1.116854 0.8101214 0.694529 1.0466708 0.3551375
    ##          [,8]      [,9]     [,10]     [,11] [,12]
    ## [1,] 1.087794 0.8277188 1.0779756 0.5651250     1
    ## [2,] 1.072401 0.7409449 0.4494975 2.8441674     1
    ## [3,] 4.888614 0.2226610 6.6084202 4.6273951     1
    ## [4,] 1.060282 0.1368132 0.1094236 1.1665615     1
    ## [5,] 1.644214 2.2686004 2.8051097 1.9136866     1
    ## [6,] 1.056410 0.3154483 0.7324739 0.5862605     1

``` r
.info <- NULL
if (is.null(.info)) {
    info_init <- list()
} else {
    info_init <- .info
  }

 deflt <- list(
                    inita    = pracma::rand(500, 500),
                    maxitr   = 30000,
                    printitr = 1e12,
                    conv     = 1e-10
                  )
  
  .info  <- initstruct(info_init, deflt)
  .info
      # for (i  in 1) {
      #   pracma::fprintf('Iter:   RMS:   Wt RMS:   Pigment Ratio RMS:   Weighted PR RMS:     dRMS:   dTaxa RMS:   dPig RMS:\n')
      #   pracma::fprintf('%5i%#7.3g%#10.3g%#21.3g%#19.3g%#10.3e%13.3e%#12.3e\n', 
      #                   itr,rms_pig,rmsxwt,rms_pig_r,rmsbwt,rms_chg,taxa_amt_chg,pig_r_chg)        # pracma::fprintf("%14 i%10.3g", 105, rms_pig)
      #   # pracma::fprintf("%14i", 2)
      # }
      # 

# testa <- info$inita  
# testb <- info$initb
# 
# dim(testa)
# dim(testb)
# 
# s - (testa %*%  testb)
```
