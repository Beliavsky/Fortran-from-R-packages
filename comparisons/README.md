# R versus Fortran comparisons

This directory contains deterministic correctness and speed comparisons between
installed R packages and their translated Fortran packages in this repository.
The shared `common` FPM package reads dated asset-price panels and constructs
explicit simple or logarithmic return panels outside timed sections. Its date
type also supports calendar fields needed by seasonal comparisons.
The suite contains 62 cases spanning probability distributions, descriptive
statistics, numerical analysis, fractional time-series operations, clustering,
geometry, polynomials, and signal processing. The existing
`mixtools/run_mvnormal_comparison.*` adds an iterative mixture-model comparison.

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
