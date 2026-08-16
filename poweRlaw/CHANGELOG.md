# Changelog

## 0.1.0

- Initial modern Fortran/FPM computational translation of `poweRlaw` 1.0.0.
- Added all eight supported fitted distribution families.
- Added direct continuous/discrete power-law d/p/r routines.
- Added parameter MLE and candidate-grid evaluation.
- Added KS/reweighted `xmin` estimation.
- Added Vuong distribution comparison.
- Added ordinary and semi-parametric bootstrap workflows and simulation helpers.
- Vendored the supplied `pracma-fortran` dependency.
- Added fast Euler-Maclaurin Hurwitz-zeta evaluation for discrete power-law fits.
- Added strict runtime-checked and optimized regression tests.
- Omitted R plotting/reference-class/presentation infrastructure.
