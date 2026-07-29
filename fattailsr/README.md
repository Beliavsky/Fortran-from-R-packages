# FatTailsR modern Fortran

A modern Fortran/FPM translation of the computational core of **FatTailsR
2.0.1**, Patrice Kiener's package for Kiener distributions and fat-tail
analysis.

The project preserves the original **GPL-2.0-only** license. The complete
license is in `COPYING`; original package metadata and citation information are
under `original/`.

## Implemented numerical scope

- Standardized logistic density, CDF, quantile, random generation, derivatives,
  left/right tail means, VaR, and expected shortfall.
- Kiener K1, K2, K3, K4, and canonical seven-parameter representation.
- Density, CDF, quantile, logit-domain quantile/inverse, probability and
  quantile derivatives, and random generation.
- Left-tail mean, right-tail mean, expected shortfall, VaR, tail-mean-minus-
  quantile, and the `c` and `h` tail correction coefficients.
- All algebraic conversions among the `a`, `k`, `w`, `d`, and `e` tail
  parameterizations used by FatTailsR.
- Raw and central moments, mean, variance, standard deviation, skewness,
  kurtosis, and excess kurtosis.
- Five-, seven-, and eleven-quantile parameter estimators.
- A bounded K4 quantile-regression fitter for vectors of observations.
- `kashp`, `ashp`, their derivative, beta/incomplete-beta support, normal
  density/CDF/quantile utilities, elevation of negative prices, missing-value
  carry-forward, and log/simple price returns.

## Deliberately omitted R infrastructure

The translation does not reproduce plotting, S3/S4 classes, data-frame/list/
`xts`/`zoo`/`timeSeries` dispatch, parallel R workers, bundled `.rda` datasets,
R help-page machinery, or website/data-download helpers. Fortran arrays and
explicit derived types replace R's dynamic containers.

The R implementation uses `minpack.lm::nlsLM` in its regression routines. The
Fortran `fit_kiener_k4` routine uses a bounded derivative-free pattern search,
initialized by the translated five-quantile estimator. This avoids an external
nonlinear least-squares dependency while fitting the same K4 quantile curve.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example fit_example
```

The package uses only standard Fortran and has no external library dependency.

## Minimal use

```fortran
program example
   use fattailsr, only : dp, kiener_parameters, make_k4, qkiener, eskiener
   implicit none
   type(kiener_parameters) :: par

   par = make_k4(m=0.0_dp, g=1.0_dp, k=6.0_dp, e=0.1_dp)
   print *, qkiener(0.01_dp, par)
   print *, eskiener(0.01_dp, par)
end program example
```

## Main modules

- `fattailsr_math`: numerical primitives and special functions.
- `fattailsr_params`: parameter type and conversion formulas.
- `fattailsr_distributions`: logistic and Kiener distributions and tail risk.
- `fattailsr_moments`: theoretical and sample moments.
- `fattailsr_estimation`: quantile estimators and K4 fitting.
- `fattailsr_returns`: price transformations and returns.
- `fattailsr`: umbrella module exporting the public API.

## Numerical conventions

- Double precision is `dp = kind(1.0d0)`.
- The version 2.x FatTailsR standardization is used, with scale factor
  `sqrt(3)/pi`, so the logistic limit has standard deviation `g`.
- Tail means require valid beta-function parameters. In practical terms, a
  corresponding tail exponent must exceed 1. Higher theoretical moments exist
  only when their order is smaller than both `a` and `w`.
- Probabilities passed to logit-based functions are clamped by machine epsilon
  away from zero and one.

## Attribution

Original package:

- Patrice Kiener, *FatTailsR: Kiener Distributions and Fat Tails in Finance and
  Neuroscience*, version 2.0.1, 2026.
- Original citation metadata is retained in `original/CITATION`.

This translation is not an official release by the original author.
