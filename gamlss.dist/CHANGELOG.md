# Changelog

## 0.3.0 - 2026-08-14

- Added ST3C, SN1, SN2, standardized skew-t (SST), generalized-t (GT), and ex-Gaussian families.
- Added GAMLSS Pareto, Pareto-I, Pareto-II, and Pareto-II-original parameterizations.
- Added double-binomial (DBI) with exact finite-support normalization.
- Added PIG2.
- Added zero-inflated and zero-adjusted PIG and Sichel families.
- Added zero-inflated/adjusted beta-binomial and beta-negative-binomial families.
- Added zero-adjusted Zipf.
- Added flexible gamma (`GAF`), flexible negative binomial (`NBF`), and zero-inflated NBF (`ZINBF`).
- Extended `fit_gamlss` from 41 to 62 supported family constants.
- Added fixed-denominator `fit_dbi`, `fit_zibb`, and `fit_zabb` likelihood fitters.
- Added v0.3 scalar-reference, normalization, round-trip, and regression-recovery tests.
- Added `example/v03_remaining.f90`.
- Retained the selected upstream R/ST3 C sources used as v0.3 algorithm references.
- Kept all Fortran source/test/example lines within the standard 132-column free-form limit.

## 0.2.0 - 2026-08-14

- Added GIG with arbitrary-real-order modified Bessel K numerical support.
- Added SHASHo/SHASH and SIMPLEX.
- Added the complete SEP through SEP4 skew-power-exponential catalog.
- Added the complete ST1 through ST5 skew-Student-t catalog.
- Added NET.
- Added generalized Poisson and double Poisson.
- Added Delaporte, SI, and Sichel count distributions using upstream recurrences.
- Added Yule, Waring, and Zipf distributions.
- Extended generic `fit_gamlss` from 18 to 41 supported family constants.
- Added independent numerical-reference tests, normalization tests, discrete
  support tests, CDF/quantile round trips, and GPO regression recovery.
- Added `example/v02_extended.f90`.
- Retained selected upstream R sources used as v0.2 algorithm references.
- Kept all Fortran source/test/example lines within the standard 132-column
  free-form limit.

## 0.1.0 - 2026-08-14

- Initial modern Fortran/FPM computational translation of gamlss.dist 6.1-1.
- Added roughly four dozen GAMLSS family parameterizations with d/p/q/r APIs.
- Added BCCG, BCT, and BCPE positive-support Box-Cox families.
- Ported the upstream Poisson-inverse-Gaussian recurrence.
- Added common GAMLSS links.
- Added generic multi-parameter maximum-likelihood regression fitting for 18
  translated families, with numerical-Hessian covariance matrices and AIC.
- Added core numerical tests and example.
- Preserved upstream metadata, reference C sources, and GPL licensing material.
