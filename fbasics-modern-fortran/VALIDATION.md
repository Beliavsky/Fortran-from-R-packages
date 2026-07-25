# Validation

## Environment

- Compiler: GNU Fortran 14.2.0
- Linear algebra: system LAPACK and BLAS
- Platform: Linux x86-64

## Debug flags

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-Wno-unused-function -Wno-compare-reals
-fcheck=all -fbacktrace -ffree-line-length-none
```

## Optimized flags

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror
-Wno-unused-function -Wno-compare-reals
-fbacktrace -ffree-line-length-none
```

## Validation command

```sh
make check
```

## Test suites

### Matrix, lag, RNG, and statistics

Tests LAPACK inversion, SVD rank, norms, matrix constructors, positive-definite repair, lags, polynomial distributed lags, descriptive and row statistics, L-moments, the reproducible LCG, and utility functions.

### Distribution density, moment, RNG, and fitting

Tests Normal and Student inversion; NIG and GH identities, moments, CDF/quantile paths, and simulation; standardized GH; RS-GLD; and Normal, Student, NIG, and GLD fitting.

### Interpolation, inference, and drawdown

Tests linear/bilinear/local-plane interpolation, correlation and location/scale tests, normality-test ranges, two-sample KS, Brownian expected drawdown, drawdown density/CDF, and realized drawdown.

### Extended algorithms

Tests:

- Stable Normal and Cauchy identities, stable simulation, ECF fitting, and the bounded MLE path
- FMKL logistic identities, FM5 equivalence, and all five extended GLD fitting modes
- GH/NIG robust moments and all GH-family fit wrappers
- Two-step/HAC and CUE GMM, J tests, Quadratic Spectral HAC, Andrews bandwidth, and linear restrictions
- EL, ET, and ETEL GEL and normalized implied weights
- Spline-density normalization, CDF/quantile inversion, and random generation
- Exact planar triangulation, ordinary-kriging interpolation, and unbiased kriging weights
- Shapiro-Wilk, Ansari-Bradley, Mood, Bartlett, Fligner-Killeen, Wilcoxon, Kruskal-Wallis, one-sample KS/Pearson, and adjusted Jarque-Bera paths


### GMM prewhitening and Andrews ARMA(1,1) bandwidth

Tests:

- Recovery of AR and MA coefficients from a simulated ARMA(1,1) series
- Independent reconstruction of the Andrews Quadratic Spectral ARMA(1,1) bandwidth formula
- Distinct AR(1) and ARMA(1,1) bandwidth paths
- VAR(1) coefficient recovery and exact recoloring-matrix identity
- Higher-order VAR(2) prewhitening
- Exact manual reconstruction of a prewhitened and recolored zero-lag HAC matrix
- Positive-semidefinite and symmetric prewhitened HAC output
- End-to-end two-step GMM with VAR prewhitening and automatic ARMA(1,1) bandwidth selection
- Graceful singular prewhitening fallback

## Applications executed

Both debug and optimized workflows run:

```text
demo_fbasics
distribution_example
analyze_csv example_returns.csv normal
analyze_csv example_returns.csv student
analyze_csv example_returns.csv nig
analyze_csv example_returns.csv gld
analyze_csv example_returns.csv fmkl
analyze_csv example_returns.csv fm5
analyze_csv example_returns.csv stable
analyze_csv example_returns.csv ssd
analyze_csv example_returns.csv hyp
analyze_csv example_returns.csv snig
```

## Result

```text
GPL-2.0-or-later source license checks passed.
Matrix, lag, RNG, and statistics tests passed.
Distribution density, moment, RNG, and fitting tests passed.
Interpolation, inference, and drawdown tests passed.
Extended stable, GLD, GH, GMM/GEL, spline-density, spatial, and test algorithms passed.
GMM prewhitening and Andrews ARMA(1,1) bandwidth tests passed.
```

Both strict configurations passed. `fpm.toml` is included, but `fpm` was unavailable and is not claimed as tested.
