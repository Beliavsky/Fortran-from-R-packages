# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `sgt` 2.0.
- Added SGT density, CDF, quantile, and random generation.
- Added skewed generalized-error (`q=Inf`) and uniform (`p=Inf`) limits.
- Preserved mean-centering and variance-adjustment conventions.
- Added elemental array operation for d/p/q routines.
- Added constant-parameter MLE and callback-based regression/general-model MLE.
- Added Nelder-Mead and BFGS optimization and numerical inference matrices.
- Added independent numerical, identity, tail, RNG, and MLE tests.
- Retained supplied `optimx` and `numDeriv` translations as reference-only
  snapshots because of conservative GPL-version compatibility concerns.
