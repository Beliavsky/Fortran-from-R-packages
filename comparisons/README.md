# R versus Fortran comparisons

This directory contains deterministic correctness and speed comparisons between
installed R packages and their translated Fortran packages in this repository.
The shared `common` FPM package reads dated asset-price panels and constructs
explicit simple or logarithmic return panels outside timed sections. Its date
type also supports calendar fields needed by seasonal comparisons.
The suite contains 151 cases spanning probability distributions, descriptive
statistics, numerical analysis, fractional time-series operations, fixed-parameter
GARCH filters, financial performance measures, autocorrelation and ARCH
diagnostics, covariance shrinkage, weighted quantile estimators, clustering,
geometry, polynomials, and signal processing. Thirty-five cases use the shared
asset-price fixture to check
fractional-difference, volatility, performance, and multivariate covariance
calculations on deterministic ETF returns. The existing
`mixtools/run_mvnormal_comparison.*` adds an iterative mixture-model comparison.
The 54-case `rfortran-core` section directly checks shared ordering, normal
and central Student-t/chi-square/F distributions, median/MAD, regularized
gamma/beta, transforms, quantiles, and
descriptive and rolling statistics and log-sum-exp reductions against R.

Run `comparisons\run_comparisons.bat` from the repository root on Windows, or
run `python comparisons/run_comparisons.py` on any platform.

The driver builds with `fpm --profile release`, checks explicit tolerances, and
writes `comparisons/results.csv`. Compilation and process startup are excluded
from kernel timings. Treat timings as exploratory and repeat them on an idle
machine before drawing conclusions.

The report includes overall and per-package arithmetic means, medians,
geometric means, minima, and maxima for R time, Fortran time, and the R/Fortran
ratio. It also counts Fortran wins, R wins, approximate ties (within 5%), and
measurements below timer resolution. The geometric mean of the per-case ratios
is the most useful single aggregate; raw times cover heterogeneous workloads.
