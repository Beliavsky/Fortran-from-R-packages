
<!-- README.md is generated from README.Rmd. Please edit that file -->

# frbinom

<!-- badges: start -->
<!-- badges: end -->

The goal of frbinom is to generate random variables of fractional
binomial distribution and compute its density, cumulative distribution,
and quantiles.

## Installation

You can install the development version of frbinom from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("leejeo25/frbinom")
```

## Example

10 random variables of a fractional binomial distribution.

``` r
library(frbinom)
rfrbinom(n=10, size=50, prob=.6, h=.7, c=.2)
#>  [1] 34 33 40 34 39 32 34 25 34 37
```

The probability density of the fractional binomial distribution.

``` r
dfrbinom(x=seq(0,50,1), size=50, prob=.6, h=.7, c=.2)
#>  [1] 3.309256e-02 1.431863e-03 1.557437e-03 1.698900e-03 1.858815e-03
#>  [6] 2.040225e-03 2.246763e-03 2.482766e-03 2.753423e-03 3.064954e-03
#> [11] 3.424806e-03 3.841903e-03 4.326913e-03 4.892553e-03 5.553920e-03
#> [16] 6.328817e-03 7.238051e-03 8.305650e-03 9.558921e-03 1.102823e-02
#> [21] 1.274636e-02 1.474727e-02 1.706397e-02 1.972538e-02 2.275181e-02
#> [26] 2.614900e-02 2.990078e-02 3.396057e-02 3.824255e-02 4.261379e-02
#> [31] 4.688921e-02 5.083182e-02 5.416082e-02 5.657003e-02 5.775753e-02
#> [36] 5.746564e-02 5.552719e-02 5.191025e-02 4.675122e-02 4.036456e-02
#> [41] 3.322041e-02 2.588705e-02 1.894504e-02 1.289015e-02 8.049503e-03
#> [46] 4.535348e-03 2.251572e-03 9.509789e-04 3.228462e-04 7.917274e-05
#> [51] 1.070436e-05
```

The histogram of fractional binomial random variables overlaid with its
density.

``` r
x<-rfrbinom(10000, 50, .6, .7, .2)
den.<-dfrbinom(seq(0,50,1), 50, .6, .7, .2)
hist(x, breaks=51, probability = TRUE, main="Histogram of fractional binomial random variables and its density")
lines(seq(0,50,1), den., type="l")
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />
