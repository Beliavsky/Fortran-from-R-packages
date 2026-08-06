
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `rrcov`: Scalable Robust Estimators with High Breakdown Point

<!-- badges: start -->

[![CRAN
version](https://www.r-pkg.org/badges/version/rrcov)](https://cran.r-project.org/package=rrcov)
[![R-CMD-check](https://github.com/valentint/rrcov/workflows/R-CMD-check/badge.svg)](https://github.com/valentint/rrcov/actions)
[![downloads](https://cranlogs.r-pkg.org/badges/rrcov)](https://cran.r-project.org/package=rrcov)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/rrcov)](https://cran.r-project.org/package=rrcov)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

The package `rrcov` provides scalable robust estimators with high
breakdown point and covers a large number of robustified multivariate
analysis methods, starting with robust estimators for the multivariate
location and covariance matrix (MCD, MVE, S, MM, SD), the deterministic
versions of MCD, S and MM estimates and regularized versions (MRCD) for
high dimensions. These estimators are used to conduct robust principal
components analysis (`PcaCov()`), linear and quadratic discriminant
analysis (`Linda()`, `Qda()`), MANOVA. Projection pursuit algorithms for
PCA to be applied in high dimensions are also available (`PcaHubert()`,
`PcaGrid()` and `PcaProj()`).

## Installation

The `rrcov` package is on CRAN (The Comprehensive R Archive Network) and
the latest release can be easily installed using the command

    install.packages("rrcov")
    library(rrcov)

## Building from source

To install the latest stable development version from GitHub, you can
pull this repository and install it using

    ## install.packages("remotes")
    remotes::install_github("valentint/rrcov", build_opts = c("--no-build-vignettes"))

Of course, if you have already installed `remotes`, you can skip the
first line (I have commented it out).

## Example

This is a basic example which shows you if the package is properly
installed:

``` r

library(rrcov)
#> Loading required package: robustbase
#> Scalable Robust Estimators with High Breakdown Point (version 1.7-3)
data(hbk)
(out <- CovMcd(hbk))
#> 
#> Call:
#> CovMcd(x = hbk)
#> -> Method:  Fast MCD(alpha=0.5 ==> h=40); nsamp = 500; (n,k)mini = (300,5) 
#> 
#> Robust Estimate of Location: 
#>       X1        X2        X3         Y  
#>  1.55833   1.80333   1.66000  -0.08667  
#> 
#> Robust Estimate of Covariance: 
#>     X1        X2        X3        Y       
#> X1   1.58739   0.03129   0.21694   0.10748
#> X2   0.03129   1.60733   0.25612   0.02864
#> X3   0.21694   0.25612   1.47254  -0.18174
#> Y    0.10748   0.02864  -0.18174   0.44081
```

## Community guidelines

### Report issues and request features

If you experience any bugs or issues or if you have any suggestions for
additional features, please submit an issue via the
[*Issues*](https://github.com/valentint/rrcov/issues) tab of this
repository. Please have a look at existing issues first to see if your
problem or feature request has already been discussed.

### Contribute to the package

If you want to contribute to the package, you can fork this repository
and create a pull request after implementing the desired functionality.

### Ask for help

If you need help using the package, or if you are interested in
collaborations related to this project, please get in touch with the
package maintainer.
