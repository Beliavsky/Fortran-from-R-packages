
[![Build Status](https://app.travis-ci.com/AnthonyChristidis/RPESE.svg?branch=master)](https://app.travis-ci.com/AnthonyChristidis/RPESE) 
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RPESE)](https://cran.r-project.org/package=RPESE)
[![Downloads](https://cranlogs.r-pkg.org/badges/RPESE)](https://cran.r-project.org/package=RPESE)

RPESE
=====

This package provides functions for computing standared error estimates for risk and performance measures of asset or portfolio returns.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=RPESE).

``` r
install.packages("RPESE", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/RPESE).

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/RPESE")
```

### Usage

``` r
# Sample Code
library(RPESE)
# Loading hedge fund data
data(edhec, package = "PerformanceAnalytics")
colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
                    "GM", "LS", "MA", "RV", "SS", "FoF")
# Computing the standard errors for the three influence functions based approaches
ES.out <- ES.SE(edhec, se.method = c("IFiid","IFcor","IFcorPW"),
                cleanOutliers = TRUE, 
                fitting.method = c("Exponential", "Gamma")[1])
# Print output
printSE(ES.out)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
