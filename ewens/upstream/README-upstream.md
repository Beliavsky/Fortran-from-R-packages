
# ewens

<!-- badges: start -->
<!-- badges: end -->

The `R` package `ewens` gives functions for the [Ewens sampling distribution](https://en.wikipedia.org/wiki/Ewens%27s_sampling_formula). The Ewens distribution is a probability distribution over partitions of an integer $n$. It is governed by a single non-negative parameter, $\theta$. The larger the value of $\theta$, the more parts there are in the partition. If $\theta$ is equal to zero, there is just a single partition with size $n$. 

The Ewens distribution can be thought of as a distribution over partitions of an integer, but it originated in population genetics to give the probability of alleles in a sample, and sampling from the Ewens distribution involves running through the [Chinese Restaurant Process](https://en.wikipedia.org/wiki/Chinese_restaurant_process) $n$ times. For further information on the Ewens sampling formula, see Simon Tavaré's [The magical Ewens sampling formula](https://www.columbia.edu/cu/simontavare/STpapers-pdf/t20b.pdf) or Harry Crane's [The Ubiquitous Ewens Sampling Formula](https://www.jstor.org/stable/24780825).

## Installation

You can install the development version of `ewens` like so:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("chrishanretty/ewens")
```

## Drawing from the Ewens distribution

Suppose I wish to divide a class of 24 into groups of different sizes. I set $\theta$ to one, and call `rewens`

``` r
set.seed(2923)
library(ewens)

allocation <- rewens(24, theta = 1)
```

The vector `allocation` gives the groups to which each student has been assigned. The sizes of the groups can be recovered by tabulating the result.

``` r
table(allocation)
```

## The probability mass function of the Ewens distribution

Suppose a student challenges their allocation into a particular group. What can we say about the probability of the particular allocation, given $\theta = 1$?

``` r
pr <- dewens(allocation, theta = 1)
pr
```

The probability of drawing this particular allocation is very low --
but then again, given that there are 1575 different partitions of
twenty-four, we should not be too surprised.

## The number of parts in the partition

It is also possible to calculate the probability that there are $x$
parts in a partition of $n$ given $\theta$. This is done using
`dewens_k`. Thus, the probability of creating the groups above is
given by

```{r}
dewens_k(k = length(unique(allocation)),
         n = 24,
         theta = 1)
```


