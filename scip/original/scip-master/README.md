
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SCIP <img src="man/figures/logo.png" width="120" align="right" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/scip)](https://CRAN.R-project.org/package=scip)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/scip)](https://CRAN.R-project.org/package=scip)
[![R-CMD-check](https://github.com/bnaras/scip/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bnaras/scip/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is an R interface to the [SCIP Optimization
Suite](https://arxiv.org/abs/2511.18580).
[SCIP](https://www.scipopt.org/) is one of the fastest non-commercial
solvers for mixed integer programming (MIP) and mixed integer nonlinear
programming (MINLP). It is also a framework for constraint integer
programming and branch-cut-and-price. It allows for total control of the
solution process and the access of detailed information down to the guts
of the solver.

## Installation

Install the released version from CRAN:

``` r
install.packages("scip")
```

Or install the development version from GitHub:

``` r
# install.packages("pak")
pak::pak("bnaras/scip")
```

## Usage

The easiest way to use SCIP (and many other solvers) is via
[`CVXR`](https://cvxr.rbind.io/) (version 1.8.2 and higher). However,
you are also welcome to refer to the [package
vignette](https://bnaras.github.io/scip/articles/scip-examples.html) for
examples.

## License

Apache License 2.0
