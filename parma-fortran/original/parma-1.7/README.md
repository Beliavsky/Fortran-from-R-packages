
[![R-CMD-check](https://github.com/alexiosg/parma/workflows/R-CMD-check/badge.svg)](https://github.com/alexiosg/parma/actions)
[![Last-changedate](https://img.shields.io/badge/last%20change-2022--10--27-yellowgreen.svg)](/commits/master)
[![packageversion](https://img.shields.io/badge/Package%20version-1.7-orange.svg?style=flat-square)](commits/master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/parma)](https://cran.r-project.org/package=parma)

# parma

The **p**ortfolio **a**llocation and **r**isk **m**anagament
**a**pplications (**parma**) package contains a unique set of methods
and models for the optimal allocation of capital in financial
portfolios. It uniquely represents certain discontinuous problems using
their smooth approximation counterparts and implements fractional based
programming for the direct optimization of risk-to-reward ratios. In
combination with the [rmgarch](https://github.com/alexiosg/rmgarch)
package, it enables the confident solution to scenario based
optimization problems using such risk and deviation measures as Mean
Absolute Deviation (MAD), Variance (EV), Minimax, Conditional Value at
Risk (CVaR), Conditional Drawdown at Risk (CDaR) and Lower Partial
Moments (LPM). In addition, it implements moment based optimization for
use with the quadratic EV problem, and a higher moment CARA utility
expansion using the coskewness and cokurtosis matrices generated from
the GO-GARCH with affine GH or NIG distributions. Benchmark relative
optimization (tracking error) is also implemented as are basic mixed
integer cardinality constraints. Finally, for non-convex problem
formulations such as the upper to lower partial moments function, global
optimization methods using a penalty based method are available.

The stable version is on
[CRAN](https://cran.r-project.org/web/packages/parma/index.html).
