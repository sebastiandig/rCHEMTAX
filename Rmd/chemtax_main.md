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

![
a = \\frac {a \\times(x \\cdot b^T)}
          {a \\cdot(b \\cdot b^T) + 10^{-100}} 
\\ ; \\ 
b = \\frac{b \\times (a^T \\cdot x + b_0)}
         {(a^T \\cdot(a \\cdot b) + b) + 10^{-100}}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Aa%20%3D%20%5Cfrac%20%7Ba%20%5Ctimes%28x%20%5Ccdot%20b%5ET%29%7D%0A%20%20%20%20%20%20%20%20%20%20%7Ba%20%5Ctimes%28b%20%5Ccdot%20b%5ET%29%20%2B%2010%5E%7B-100%7D%7D%20%0A%5C%20%3B%20%5C%20%0Ab%20%3D%20%5Cfrac%7Bb%20%5Ctimes%20%28a%5ET%20%5Ccdot%20x%20%2B%20b_0%29%7D%0A%20%20%20%20%20%20%20%20%20%7B%28a%5ET%20%5Ccdot%28a%20%5Ccdot%20b%29%20%2B%20b%29%20%2B%2010%5E%7B-100%7D%7D%0A "
a = \frac {a \times(x \cdot b^T)}
          {a \times(b \cdot b^T) + 10^{-100}} 
\ ; \ 
b = \frac{b \times (a^T \cdot x + b_0)}
         {(a^T \cdot(a \cdot b) + b) + 10^{-100}}
")

Iterations with standard deviations:

![
\\begin{aligned}
& a\_{new} = \\frac {a \\times(x\_{w} \\cdot b^T)}
                {((a \\cdot b ) \\times x\_{\\sigma^2}^{-1}) \\cdot b^T + 10^{-100}} 
\\ ; \\ 
b\_{new} = \\frac {b \\times (a^T \\cdot x_w) + b\_{0w}} 
                {a^T \\cdot ((a \\cdot b) \\times x\_{\\sigma^2}^{-1}) + b \\times b\_{\\sigma^2}^{-1} + 10^{-100}} \\\\\\\\
& \\text {where}; \\\\
& x = \\text{pigment concentations for each sample}; \\\\
& x\_{\\sigma^2}^{-1} = \\frac {1}{x\_{\\sigma}^2} = \\text{inverse-variance weighting of }\\ x; \\\\ 
& x_w = x \\times x\_{\\sigma^2}^{-1} = \\text{weighted } x \\text{; } \\\\
& a = \\text{taxa contribution}; \\\\
& b = \\text{pigment ratio matrix};\\\\ 
& b\_{\\sigma^2}^{-1} = \\frac {1}{b\_{0\\sigma}^2} = \\text{inverse-variance weighting of initial}\\ b; \\\\ 
& b\_{0w} = b_0 \\times b\_{\\sigma^2}^{-1}  = \\text{weighted initial}\\ b \\\\
\\end{aligned}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Baligned%7D%0A%26%20a_%7Bnew%7D%20%3D%20%5Cfrac%20%7Ba%20%5Ctimes%28x_%7Bw%7D%20%5Ccdot%20b%5ET%29%7D%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7B%28%28a%20%5Ccdot%20b%20%29%20%5Ctimes%20x_%7B%5Csigma%5E2%7D%5E%7B-1%7D%29%20%5Ccdot%20b%5ET%20%2B%2010%5E%7B-100%7D%7D%20%0A%5C%20%3B%20%5C%20%0Ab_%7Bnew%7D%20%3D%20%5Cfrac%20%7Bb%20%5Ctimes%20%28a%5ET%20%5Ccdot%20x_w%29%20%2B%20b_%7B0w%7D%7D%20%0A%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%20%7Ba%5ET%20%5Ccdot%20%28%28a%20%5Ccdot%20b%29%20%5Ctimes%20x_%7B%5Csigma%5E2%7D%5E%7B-1%7D%29%20%2B%20b%20%5Ctimes%20b_%7B%5Csigma%5E2%7D%5E%7B-1%7D%20%2B%2010%5E%7B-100%7D%7D%20%5C%5C%5C%5C%0A%26%20%5Ctext%20%7Bwhere%7D%3B%20%5C%5C%0A%26%20x%20%3D%20%5Ctext%7Bpigment%20concentations%20for%20each%20sample%7D%3B%20%5C%5C%0A%26%20x_%7B%5Csigma%5E2%7D%5E%7B-1%7D%20%3D%20%5Cfrac%20%7B1%7D%7Bx_%7B%5Csigma%7D%5E2%7D%20%3D%20%5Ctext%7Binverse-variance%20weighting%20of%20%7D%5C%20x%3B%20%5C%5C%20%0A%26%20x_w%20%3D%20x%20%5Ctimes%20x_%7B%5Csigma%5E2%7D%5E%7B-1%7D%20%3D%20%5Ctext%7Bweighted%20%7D%20x%20%5Ctext%7B%3B%20%7D%20%5C%5C%0A%26%20a%20%3D%20%5Ctext%7Btaxa%20contribution%7D%3B%20%5C%5C%0A%26%20b%20%3D%20%5Ctext%7Bpigment%20ratio%20matrix%7D%3B%5C%5C%20%0A%26%20b_%7B%5Csigma%5E2%7D%5E%7B-1%7D%20%3D%20%5Cfrac%20%7B1%7D%7Bb_%7B0%5Csigma%7D%5E2%7D%20%3D%20%5Ctext%7Binverse-variance%20weighting%20of%20initial%7D%5C%20b%3B%20%5C%5C%20%0A%26%20b_%7B0w%7D%20%3D%20b_0%20%5Ctimes%20b_%7B%5Csigma%5E2%7D%5E%7B-1%7D%20%20%3D%20%5Ctext%7Bweighted%20initial%7D%5C%20b%20%5C%5C%0A%5Cend%7Baligned%7D%0A "
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
")

Non-negative Least Squares:

-   Minimize
    ![a](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;a "a")
    matrix

![
\\begin{aligned}
& min\|\|x_i - \\hat{a_i} \\cdot b_i\|\| \\\\
& i = \\text{sample row}; \\\\
& x_i = (\\frac {x_i} {x\_{\\sigma_i}})^2; \\\\
& b_i = \\frac {b^T} {x\_{\\sigma_i}}; \\\\
& \\hat{a_i} = \\text{estimated taxa contribution per sample}\\\\
\\end{aligned}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Baligned%7D%0A%26%20min%7C%7Cx_i%20-%20%5Chat%7Ba_i%7D%20%5Ccdot%20b_i%7C%7C%20%5C%5C%0A%26%20i%20%3D%20%5Ctext%7Bsample%20row%7D%3B%20%5C%5C%0A%26%20x_i%20%3D%20%28%5Cfrac%20%7Bx_i%7D%20%7Bx_%7B%5Csigma_i%7D%7D%29%5E2%3B%20%5C%5C%0A%26%20b_i%20%3D%20%5Cfrac%20%7Bb%5ET%7D%20%7Bx_%7B%5Csigma_i%7D%7D%3B%20%5C%5C%0A%26%20%5Chat%7Ba_i%7D%20%3D%20%5Ctext%7Bestimated%20taxa%20contribution%20per%20sample%7D%5C%5C%0A%5Cend%7Baligned%7D%0A "
\begin{aligned}
& min||x_i - \hat{a_i} \cdot b_i|| \\
& i = \text{sample row}; \\
& x_i = (\frac {x_i} {x_{\sigma_i}})^2; \\
& b_i = \frac {b^T} {x_{\sigma_i}}; \\
& \hat{a_i} = \text{estimated taxa contribution per sample}\\
\end{aligned}
")

-   Minimize
    ![b](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;b "b")
    matrix

    ![
    \\begin{aligned}
    & min\|\|x\_{ij} - a\_{ik} \\cdot \\hat{b\_{jk}}\|\| \\\\
    & i = \\text{sample};\\
    j = \\text{pigment};\\
    k = \\text{taxa}; \\\\
    & k,j \\in \\Phi; \\ 
    \\Phi = \\text{set of}\\ k,j\\ \\text{non zero elements in pigment ratio matrix} \\\\
    & x\_{j} = 
    \\begin{pmatrix} 
            \\frac {x\_{ij}} {x\_{\\sigma\_{ij}}} \\\\
            \\frac {b\_{jk}} {b\_{\\sigma\_{jk}}}
    \\end{pmatrix} ;\\\\
    & a\_{kj} = 
    \\begin{pmatrix}  
            \\frac {a\_{ijk}} {x\_{\\sigma\_{ij}}} \\\\
            diag( b\_{\\sigma\_{kj}}^{-1} )
    \\end{pmatrix} ;\\\\
    & \\hat {b\_{jk}} = \\text {estimated pigment ratio matrix} \\\\
    \\end{aligned}
    ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cbegin%7Baligned%7D%0A%26%20min%7C%7Cx_%7Bij%7D%20-%20a_%7Bik%7D%20%5Ccdot%20%5Chat%7Bb_%7Bjk%7D%7D%7C%7C%20%5C%5C%0A%26%20i%20%3D%20%5Ctext%7Bsample%7D%3B%5C%0Aj%20%3D%20%5Ctext%7Bpigment%7D%3B%5C%0Ak%20%3D%20%5Ctext%7Btaxa%7D%3B%20%5C%5C%0A%26%20k%2Cj%20%5Cin%20%5CPhi%3B%20%5C%20%0A%5CPhi%20%3D%20%5Ctext%7Bset%20of%7D%5C%20k%2Cj%5C%20%5Ctext%7Bnon%20zero%20elements%20in%20pigment%20ratio%20matrix%7D%20%5C%5C%0A%26%20x_%7Bj%7D%20%3D%20%0A%5Cbegin%7Bpmatrix%7D%20%0A%20%20%20%20%20%20%20%20%5Cfrac%20%7Bx_%7Bij%7D%7D%20%7Bx_%7B%5Csigma_%7Bij%7D%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20%5Cfrac%20%7Bb_%7Bjk%7D%7D%20%7Bb_%7B%5Csigma_%7Bjk%7D%7D%7D%0A%5Cend%7Bpmatrix%7D%20%3B%5C%5C%0A%26%20a_%7Bkj%7D%20%3D%20%0A%5Cbegin%7Bpmatrix%7D%20%20%0A%20%20%20%20%20%20%20%20%5Cfrac%20%7Ba_%7Bijk%7D%7D%20%7Bx_%7B%5Csigma_%7Bij%7D%7D%7D%20%5C%5C%0A%20%20%20%20%20%20%20%20diag%28%20b_%7B%5Csigma_%7Bkj%7D%7D%5E%7B-1%7D%20%29%0A%5Cend%7Bpmatrix%7D%20%3B%5C%5C%0A%26%20%5Chat%20%7Bb_%7Bjk%7D%7D%20%3D%20%5Ctext%20%7Bestimated%20pigment%20ratio%20matrix%7D%20%5C%5C%0A%5Cend%7Baligned%7D%0A "
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
    ")

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

-   `length(x)` = largest of the dimensions, so if 100 x 4 = 100, or 4 x
    100 = 100
-   `*` is matrix multiply,
-   `.*` & `./` are matrix element by element multiply & divide
    operations
-   can use `pracma::rand()` for rand
-   `repmat` is used to replicate a matrix by `x * y`, ie
    `repmat({matrix 5 x 4}, 2, 4)` result would be matrix 2 x 16, so 8
    replicates of the original matrix where stacked twice and multiplied
    by 4
-   `mod(a,m)` is same as `a %% m` in R
-   `df(:)` puts matrix into 1 column
-   `df(end)` gets the last digit, bottom right of matrix

## Normalization to Chlor-a

I think this is done in `matfactuvw` using `pracma::Diag(1/chlor-a col)`
then the dot product with df

\#TODO: something to include is: This would load internal ratios file so
that you have a basic pigment ratio var
`pigment_ratio <- system.file("", package = "rCHEMTAX")`

## Plot in log notation on axes

\# allow y or x log style plotting from Matlab in semilogy, semilogx,
loglog

![
\\text{i.e. }x\\times10^{n}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Ctext%7Bi.e.%20%7Dx%5Ctimes10%5E%7Bn%7D%0A "
\text{i.e. }x\times10^{n}
")

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

-   `repmat`
-   `fprintf`
-   `ceil` - not anymore, was alias for `base::ceiling`
-   `rand`
-   `Diag` - not anymore
-   `lsqnonneg` - could replace with `limsolve`

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

    ##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
    ## [1,] 0.2571389 0.9806986 0.7195641 0.1371409 0.7144976 0.7137136 0.6168906
    ## [2,] 0.2229049 0.1436677 1.5336944 0.1106050 0.6688343 1.8588924 1.7933113
    ## [3,] 0.5855983 2.1998345 2.3833450 1.4809635 0.3884695 2.4787525 1.5631827
    ## [4,] 0.7439398 1.4553259 1.6961080 0.1891135 0.4789913 0.7586944 0.3695204
    ## [5,] 0.5017947 1.7463629 4.8506696 2.3446969 3.4920936 1.8773109 4.2311748
    ## [6,] 1.9021084 2.6892627 3.4833697 1.5772124 1.8363765 2.7513222 1.4575227
    ##           [,8]       [,9]     [,10]      [,11] [,12]
    ## [1,] 0.3555101 0.02602217 0.5337614 1.22371100     1
    ## [2,] 1.0942676 0.27885945 0.5369586 0.11412349     1
    ## [3,] 0.1362293 1.01472824 2.2746197 0.92638869     1
    ## [4,] 0.3850744 1.21898753 1.4966944 0.05756243     1
    ## [5,] 2.0081428 2.04126945 3.0444069 0.03012610     1
    ## [6,] 3.2899316 2.83728970 0.7367987 1.54661710     1

``` r
# apply method *more confusing
t(apply(test_norm, 1, function(x) x / x[length(x)]))
```

    ##           [,1]      [,2]      [,3]      [,4]      [,5]      [,6]      [,7]
    ## [1,] 0.2571389 0.9806986 0.7195641 0.1371409 0.7144976 0.7137136 0.6168906
    ## [2,] 0.2229049 0.1436677 1.5336944 0.1106050 0.6688343 1.8588924 1.7933113
    ## [3,] 0.5855983 2.1998345 2.3833450 1.4809635 0.3884695 2.4787525 1.5631827
    ## [4,] 0.7439398 1.4553259 1.6961080 0.1891135 0.4789913 0.7586944 0.3695204
    ## [5,] 0.5017947 1.7463629 4.8506696 2.3446969 3.4920936 1.8773109 4.2311748
    ## [6,] 1.9021084 2.6892627 3.4833697 1.5772124 1.8363765 2.7513222 1.4575227
    ##           [,8]       [,9]     [,10]      [,11] [,12]
    ## [1,] 0.3555101 0.02602217 0.5337614 1.22371100     1
    ## [2,] 1.0942676 0.27885945 0.5369586 0.11412349     1
    ## [3,] 0.1362293 1.01472824 2.2746197 0.92638869     1
    ## [4,] 0.3850744 1.21898753 1.4966944 0.05756243     1
    ## [5,] 2.0081428 2.04126945 3.0444069 0.03012610     1
    ## [6,] 3.2899316 2.83728970 0.7367987 1.54661710     1

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

$$

$$
