# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational port of REBayes 2.60.
- Replaced Rmosek/MOSEK finite-grid KW optimization with native EM plus
  vertex-direction acceleration and KKT diagnostics.
- Added principal one- and two-dimensional mixture families, repeated-measures
  mixtures, posterior/FDR utilities, RLR, MEDDE and numerical utilities.
- Added deterministic noncentral-t quadrature.
- Added strict runtime/interface test suite and independent randomized
  validation.
- Plotting, R object infrastructure, Rxiv, BDGLmix, HuberSpline and
  HodgesLehmann are not included in v0.1.0.
