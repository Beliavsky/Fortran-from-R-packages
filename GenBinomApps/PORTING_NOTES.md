# Porting notes

## Scope

This is a complete translation of the computational code in `GenBinomApps`
1.2.1. The package contains no plotting or compiled native-code layer.

## Generalized binomial distribution

The upstream implementation expands every grouped binomial into Bernoulli
trials and repeatedly convolves two-point distributions. The Fortran routine
uses the same exact Bernoulli-convolution recurrence in-place. A PMF or CDF
table is built once for vector evaluation, quantiles, and random generation.

The distribution is also known as the Poisson-binomial distribution.

## Beta functions

Base R's `pbeta` and `qbeta` are replaced by a standalone regularized
incomplete-beta implementation (continued fractions) and inverse-beta
bisection. This is sufficient for the exact Clopper-Pearson formulas used by
the package and avoids an external special-function dependency.

The countermeasure lower-limit formula in the R source uses
`pbeta(p,1e-100,1+n)` as a proxy for the limiting Beta(0,b) CDF. The Fortran
implementation evaluates that intended limiting behavior directly: zero at
`p=0`, one at every positive `p`.

## Sample-size searches

The R routines solve a continuous equation with `uniroot` and return
`ceiling(root)`, with a default artificial upper bound of `1e100`.

The Fortran routines instead bracket and binary-search the smallest **integer**
sample size satisfying the same monotone Clopper-Pearson inequality. This is
the integer returned by `ceiling(root)` but is more robust and avoids
evaluating beta functions with meaningless parameters near `1e100`.

An optional `max_n` integer bound replaces R's `uniroot.upper`.

## RNG

`rgbinom` retains inversion sampling from the exact generalized-binomial CDF.
Fortran's intrinsic RNG is used, so seeded streams do not reproduce R's
`runif()` bit-for-bit.
