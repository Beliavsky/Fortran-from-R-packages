# BivKLD

`BivKLD` estimates directed Kullback-Leibler (KL) divergence between
bivariate distributions. Its kernel estimator follows Chackochan, Sankaran,
and Unnikrishnan Nair (2026) and uses bivariate Gaussian kernel densities from
the `ks` package.

## Installation

Install the submitted source archive with:

```r
install.packages("BivKLD_0.1.0.tar.gz", repos = NULL, type = "source")
```

## Two-sample estimate

```r
library(BivKLD)

set.seed(2026)
x <- cbind(rnorm(100), rnorm(100))
y <- cbind(rnorm(100, 0.5), rnorm(100, -0.25))

biv_kld(x, y, bandwidth = "normal")
biv_kld(y, x, bandwidth = "normal")
```

KL divergence is directed, so the two results need not be equal. The default
bandwidth selector is unconstrained smoothed cross-validation (`"scv"`), as in
the reference paper. The normal-scale selector above is useful for quick
exploration.

## Pairwise group comparisons

```r
measurements <- iris[, c("Sepal.Length", "Sepal.Width")]
biv_kld_matrix(measurements, iris$Species, bandwidth = "normal")
```

Rows are source distributions and columns are comparison distributions, so
entry `[i, j]` estimates `D(group_i || group_j)`.

## Exact calculations

```r
biv_kld_discrete(c(0.2, 0.3, 0.1, 0.4),
                 c(0.1, 0.4, 0.2, 0.3))

biv_kld_normal(c(-2, 2), matrix(c(1, 1, 1, 2), 2),
               c(-2, 2), matrix(c(3, 3, 3, 6), 2))
```

## Reference

Chackochan, R., Sankaran, P. G., & Unnikrishnan Nair, N. (2026). Bivariate
Kullback-Leibler divergence. *Communications in Statistics - Theory and
Methods*, 55(1), 292-312.
[doi:10.1080/03610926.2025.2496687](https://doi.org/10.1080/03610926.2025.2496687)

