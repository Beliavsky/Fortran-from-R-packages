# Skellam

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![CRAN](https://www.r-pkg.org/badges/version-ago/skellam)
[![R-CMD-check](https://github.com/monty-se/skellam/actions/workflows/rhub.yaml/badge.svg)](https://github.com/monty-se/skellam/actions/workflows/rhub.yaml)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/skellam)](https://cranlogs.r-pkg.org/badges/grand-total/skellam)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/skellam)](https://cran.r-project.org/package=skellam)
## Overview

The **skellam** package provides functions for working with the Skellam distribution – the distribution of the difference between two independent Poisson random variables. It includes routines for:
- Calculating the probability mass function (`dskellam`)
- Computing the cumulative distribution function (`pskellam`)
- Determining quantiles (`qskellam`)
- Generating random variates (`rskellam`)
- Performing maximum likelihood estimation (`skellam.mle`)
- Conducting regression analysis under the Skellam model (`skellam.reg`)

This package is designed to offer enhanced numerical accuracy and robust handling of a wide range of parameter values.

## Installation

Install the latest stable version from CRAN:

```r
install.packages("skellam")
```
Alternatively, install the development version from GitHub:

```r
# install.packages("remotes")  # Uncomment if needed
remotes::install_github("monty-se/skellam")
```
## Usage
### Distribution Functions

```r
dskellam(x, lambda1, lambda2 = lambda1, log = FALSE)
```
Returns the (log) density of the Skellam distribution.

```r
pskellam(q, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
```
Computes the (log) cumulative distribution function.

```r
qskellam(p, lambda1, lambda2 = lambda1, lower.tail = TRUE, log.p = FALSE)
```
Returns the quantile function for the Skellam distribution.

```r
rskellam(n, lambda1, lambda2 = lambda1)
```
Generates random variates following the Skellam distribution.

### Additional Functionalities

```r
skellam.mle(x)
```
Performs maximum likelihood estimation (MLE) for the Skellam distribution parameters based on observed differences.

```r
skellam.reg(y, x)
```
Fits a regression model assuming a Skellam distribution, using an exponential link to ensure positivity of the rate parameters.

## Theoretical Background

If $X \sim \text{Poisson}(\lambda_1)$ and $Y \sim \text{Poisson}(\lambda_2)$ are independent, then the difference

$$
Z = X - Y
$$

follows a Skellam distribution:

$$
Z \sim \text{Skellam}(\lambda_1, \lambda_2)
$$

This property is the theoretical foundation behind the functions provided in the skellam package. For more details, see the Wikipedia page on the Skellam distribution ​

## Examples

### Random Variate Generation and Density Estimation

#### Set parameters

```r
N <- 5000
lambda1 <- 1.5
lambda2 <- 0.5

```
#### Generate independent Poisson samples and compute their difference

```r
X <- rpois(N, lambda1)
Y <- rpois(N, lambda2)
XminusY <- X - Y
```
#### Generate Skellam random variates

```r
Z <- rskellam(N, lambda1, lambda2)
```

#### Plot empirical and theoretical densities

```r
xseq <- seq(min(XminusY), max(XminusY))
plot(table(XminusY), main = "Empirical vs. Theoretical Density", 
     xlab = "X - Y", ylab = "Frequency", pch = 1)
points(xseq, N * dskellam(xseq, lambda1, lambda2), col = "blue", pch = 4)
legend("topright", legend = c("Empirical", "Theoretical"), pch = c(1, 4), 
       col = c("black", "blue"))
```

### Maximum Likelihood Estimation

#### Estimate Skellam parameters from the difference data

```r
mle_result <- skellam.mle(XminusY)
print(mle_result)
```

### Skellam Regression

#### Simulate covariate data and corresponding Poisson outcomes

```r
set.seed(123)
x_cov <- rnorm(N)
y1 <- rpois(N, exp(1 + x_cov))
y2 <- rpois(N, exp(-1 + x_cov))
y_diff <- y2 - y1
```
#### Fit a Skellam regression model

```r
reg_result <- skellam.reg(y_diff, x_cov)
print(reg_result)
```
## References

Skellam, J. G. (1946). The frequency distribution of the difference between two Poisson variates belonging to different populations. Journal of the Royal Statistical Society, Series A, 109(3), 296-296.

Johnson, N. L. (1959). On an extension of the connection between Poisson and χ² distributions. Biometrika, 46, 352–362.

Wikipedia: [Skellam distribution](https://en.wikipedia.org/wiki/Skellam_distribution)

## License

The skellam package is licensed under the GPL (>= 2).

## Maintainers

Montasser Ghachem (mg@pinstimation.com)
Oguz Ersan (oe@pinstimation.com)
