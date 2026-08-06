
[![Build Status](https://app.travis-ci.com/AnthonyChristidis/RPEGLMEN.svg?branch=master)](https://app.travis-ci.com/AnthonyChristidis/RPEGLMEN) [![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/RPEGLMEN)](https://cran.r-project.org/package=RPEGLMEN) [![Downloads](https://cranlogs.r-pkg.org/badges/RPEGLMEN)](https://cran.r-project.org/package=RPEGLMEN)

RPEGLMEN
========

This package provides an implementation of the elastic net penalty for Gamma and exponentially distributed response variables.

------------------------------------------------------------------------

### Installation

You can install the **stable** version on [R CRAN](https://cran.r-project.org/package=RPEGLMEN).

``` r
install.packages("RPEGLMEN", dependencies = TRUE)
```

You can install the **development** version from [GitHub](https://github.com/AnthonyChristidis/RPEGLMEN).

``` r
library(devtools)
devtools::install_github("AnthonyChristidis/RPEGLMEN")
```

### Background Information

This package is designed to provide the user to fit an Exponential or Gamma distribution to the response variable with an elastic net penalty on the predictors. This package is of particular use in combination with the [RPEIF](https://github.com/AnthonyChristidis/RPEIF) and [RPESE](https://github.com/AnthonyChristidis/RPESE) packages, in which the influence function of a time series of returns is used to compute the standard error of a risk and performance measure. See [Chen and Martin (2018)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3085672) for more details.

For the computational details to fit a Gamma distribution on data with an elastic net penalty, see [Chen, Arakvin and Martin (2018)](https://arxiv.org/abs/1804.07780).

### Usage

``` r
# Sample Code

# Load the package
library(RPEGLMEN)

# Function to return the periodogram of data series
myperiodogram <- function (data, max.freq = 0.5, twosided = FALSE, keep = 1){
  data.fft <- fft(data)
  N <- length(data)
  tmp <- Mod(data.fft[2:floor(N/2)])^2/N
  tmp <- sapply(tmp, function(x) max(1e-05, x))
  freq <- ((1:(floor(N/2) - 1))/N)
  tmp <- tmp[1:floor(length(tmp) * keep)]
  freq <- freq[1:floor(length(freq) * keep)]
  if (twosided) {
    tmp <- c(rev(tmp), tmp)
    freq <- c(-rev(freq), freq)
  }
  return(list(spec <- tmp, freq <- freq))
}

# Function to compute the standard error based the periodogram of the influence functions time series
SE.Gamma <- function(data, d = 7, alpha = 0.5, keep = 1, exponential.dist = TRUE){
  N<-length(data)
  # Compute the periodograms
  my.periodogram <- myperiodogram(data)
  my.freq <- my.periodogram$freq
  my.periodogram <- my.periodogram$spec
  # Remove values of frequency 0 as it does not contain information about the variance
  my.freq <- my.freq[-1]
  my.periodogram <- my.periodogram[-1]
  # Implement cut-off
  nfreq <- length(my.freq)
  my.freq <- my.freq[1:floor(nfreq*keep)]
  my.periodogram <- my.periodogram[1:floor(nfreq*keep)]
  # GLM with BFGS optimization
  # Create 1, x, x^2, ..., x^d
  x.mat <- rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat <- cbind(x.mat,my.freq^col.iter)
  }
  # Fit the Exponential or Gamma model
  if(exponential.dist)
    res <- glmnet_exp(x.mat, my.periodogram, alpha.EN = alpha) else
      res <- fit.glmGammaNet(x.mat, my.periodogram, alpha.EN = alpha)
  # Return the estimated variance
  return(sqrt(exp(res[1])/N))
}

# Loading hedge fund data from PA
data(edhec, package <- "PerformanceAnalytics")
colnames(edhec)

# Computing the expected shortfall for the time series of returns
library(RPEIF)
test.mat <- apply(edhec, 2, IF.ES)
test.mat <- apply(test.mat, 2, as.numeric)

# Returning the standard errors from the Exponential distribution fit
apply(test.mat, 2, SE.Gamma, exponential.dist = TRUE)
```

### License

This package is free and open source software, licensed under GPL (&gt;= 2).
