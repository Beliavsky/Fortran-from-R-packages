# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `countDM` 0.1.0.
- Ported all exported computational routines.
- Added Touchard/Bell, Bell--Touchard and zero-inflated distribution kernels.
- Added self-contained Lambert W, BFGS optimization, numerical Hessians,
  covariance matrices, standard errors and AIC.
- Added all upstream MLE models and bundled datasets.
- Corrected the upstream `mle_zibellt` start-order defect and `qbellt`
  probability-semantics defects.
- Added strict Fortran 2018 tests and an example.
