[![Build Status](https://app.travis-ci.com/AnthonyChristidis/RPEIF.svg?branch=master)](https://app.travis-ci.com/AnthonyChristidis/RPEIF) 
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/RPEIF)](https://cran.r-project.org/package=RPEIF)
[![Downloads](https://cranlogs.r-pkg.org/badges/RPEIF)](https://cran.r-project.org/package=RPEIF)

RPEIF
=====

This package provides functions for computing the influence functions of risk and performance measures.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=RPEIF).

``` r
install.packages("RPEIF", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/RPEIF).

``` r
library(devtools) 
devtools::install_github("AnthonyChristidis/RPEIF")
```

### Usage

``` r
# Sample Code
library(RPEIF)
# Loading the data
data(edhec, package = "PerformanceAnalytics")
colnames(edhec) = c("CA", "CTAG", "DIS", "EM","EMN", "ED", "FIA",
                    "GM", "LS", "MA", "RV", "SS", "FoF")
# Computing the IF of the returns (with outlier cleaning and prewhitening)
outIF <- IF(risk = "Mean",
            returns = edhec[,"CA"], evalShape = FALSE, retVals = NULL, nuisPars = NULL,
            IFplot = TRUE, IFprint = TRUE,
            prewhiten = TRUE,
            cleanOutliers = TRUE, cleanMethod = "locScaleRob", eff = 0.99)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
