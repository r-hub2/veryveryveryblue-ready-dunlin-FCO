
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Package FCO - Flexible Cutoffs for Model Fit Evaluation in Covariance-Based Structural Models

## Changes in FCO

### FCO 0.8.0

-   Changed vignette after release on CRAN
-   Fixed issue in index_guess
-   Allows to specify sample size n instead of dataset x

------------------------------------------------------------------------

### FCO 0.7.2

-   New seed argument in gen_fit for reproducible cutoffs

### FCO 0.7.1

-   Added a `NEWS.md` file to track changes to the package.
-   Minor revisions to tests

### FCO 0.7.0

-   Speed improvements in the vignette
-   New naming scheme
-   Minor revisions to the descriptions and references
-   Added contributor

### FCO 0.69

-   Bug fixes in gen_fit for OS compatibility
-   Improvements in the vignette

### FCO 0.67

-   First stable release

## Description

The goal of FCO is to to derive flexible cutoffs for fit indices in
Covariance-based Structural Equation Modeling based on the paper by
Niemand & Mai (2018). Flexible cutoffs are an alternative to fixed
cutoffs - rules-of-thumb - regarding an appropriate cutoff for fit
indices such as CFI or SRMR. It has been demonstrated that these
flexible cutoffs perform better than fixed cutoffs in grey areas where
misspecification is not easy to detect. The package provides an
alternative to the tool at flexiblecutoffs.org as it allows to tailor
flexible cutoffs to a given dataset and model, which is so far not
available in the tool. The package simulates fit indices based on a
given dataset and model and then estimates the flexible cutoffs. Some
useful functions, e.g., to determine the GoF or BoF-nature of a fit
index, are provided. So far, additional options for a relative use (is a
model better than another?) are provided in an exploratory manner.

## Installation

You can install the FCO from CRAN [CRAN](https://cran.r-project.org)
with:

``` r
install.packages("FCO")
library(FCO)
```

## Example

This is the basic usage for FCO in case of deriving flexible cutoffs for
a single model:

``` r
library(FCO)
library(lavaan)
#> This is lavaan 0.6-11
#> lavaan is FREE software! Please report any bugs.
#Data from bb1992
mod <- "
F1 =~ Q5 + Q7 + Q8
F2 =~ Q2 + Q4
F3 =~ Q10 + Q11 + Q12 + Q13 + Q18 + Q19 + Q20 + Q21 + Q22
F4 =~ Q1 + Q17
F5 =~ Q6 + Q14 + Q15 + Q16
"

#Flexible cutoffs for this model
fits.single <- gen_fit(mod1 = mod, x = bb1992, rep = 10)
flex_co(fits = fits.single, index = c("CFI", "SRMR"))
#> Warning in flex_co(fits = fits.single, index = c("CFI", "SRMR")): The number of
#> replications is lower than the recommended minimum of 500. Consider with care.
#> $cutoff
#>        CFI       SRMR 
#> 0.97826871 0.03659316 
#> 
#> $index
#> [1] "CFI"  "SRMR"
#> 
#> $alpha
#> [1] 0.05
#> 
#> $gof
#>   CFI  SRMR 
#>  TRUE FALSE 
#> 
#> $replications
#> [1] 10
#> 
#> $`number of non-converging models`
#> [1] 0
#> 
#> $`share of non-converging models`
#> [1] 0

#Use recommend function
recommend(fits.single)
#> Warning in recommend(fits.single): The number of replications is lower than the
#> recommended minimum of 500. Consider with care.
#> $recommended
#>      type fit.values
#> SRMR  BoF      0.038
#> 
#> $cutoffs
#>               SRMR
#> cutoff 0.001 0.037
#> cutoff 0.01  0.037
#> cutoff 0.05  0.037
#> cutoff 0.1   0.036
#> 
#> $decisions
#>                  SRMR
#> cutoff 0.001 rejected
#> cutoff 0.01  rejected
#> cutoff 0.05  rejected
#> cutoff 0.1   rejected
#> 
#> $replications
#> [1] 10
#> 
#> $comment
#> [1] "Recommendations based on flexible cutoffs and Mai et al. (2021)"
```
