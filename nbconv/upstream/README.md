
# nbconv

<!-- badges: start -->
<!-- badges: end -->

## Introduction

nbconv was written to facilitate the evaluation of the sums of arbitrary negative binomial (NB) random variables. nbconv implements three distinct methods for approaching this problem: Furman's exact probability mass function (https://ssrn.com/abstract=1650365), saddlepoint approximation, and a method of moments approximation. I would like to acknowledge Martin Modrák for inspiring me to pursue the saddlepoint approximation with his [related blog post](https://www.martinmodrak.cz/2019/06/20/approximate-densities-for-sums-of-variables-negative-binomials-and-saddlepoint/).

It should be noted that Furman's PMF is expressed as a series representation. Because of this, the accuracy of the evaluation is limited by the number of terms included in the series. The probability values calculated in this way are therefore also approximate solutions, meaning that nbconv only offers approximate solutions to the convolution of NB random variables. Nevertheless, the approximations appear to be rather good ones.

Like other distribution functions in R, nbconv can calculate the density, distribution, and quantile functions of NB convolutions. These are done with <code>dnbconv()</code>, <code>pnbconv()</code>, and <code>qnbconv()</code>, respectively. Each of these functions requires specification of the method of evaluation. Random deviates can be sampled with the <code>rnbconv()</code> function and are obtained independently of any evaluation method. Finally, the function <code>nbconv_params()</code> explicitly calculates the mean, variance, skewness, and excess kurtosis of the convolution distribution using cumulants.

See the manual for more information and please do not hesitate to reach out with any questions or comments.

## Installation

You can install the development version of nbconv from [GitHub](https://github.com/gbedwell/nbconv) with:

``` r
# install.packages("devtools")
devtools::install_github("gbedwell/nbconv")
```

