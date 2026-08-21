# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of GB2 2.1.2 computational code.
- Added GB2 distribution, moments, inequality/poverty indicators, likelihood,
  analytic derivatives, Fisher information, full/profile pseudo-likelihood
  fits, empirical/robust/nonlinear-fit helpers, compound and auxiliary compound
  models, and covariance/delta-method routines.
- Integrated translated `survey-fortran` for design-based score covariance.
- Added internal adaptive Gauss-Kronrod integration and real `3F2(1)` series to
  avoid linking the supplied license-incompatible cubature/hypergeo pair.
- Added independent numerical references and derivative/fit/mixture/survey/RNG
  tests.
- Corrected documented upstream numerical/argument edge cases described in
  `PORTING_NOTES.md`.
