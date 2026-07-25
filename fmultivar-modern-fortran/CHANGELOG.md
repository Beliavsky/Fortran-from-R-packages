# Changelog

## 0.1.0

- Initial modern Fortran translation of the non-plotting computational surface
  of `fMultivar` 4031.84.
- Added bivariate Normal, Student-t, Cauchy, and seven-family elliptical density
  routines.
- Added multivariate Normal, Student-t, skew-Normal, skew-Student, and
  skew-Cauchy density, probability, quantile, and simulation workflows.
- Added Normal, skew-Normal, skew-t, and skew-Cauchy fitting with numerical
  Hessians and covariance matrices.
- Added adaptive integration, KDE, histogram, grid, square-binning, and
  hex-binning routines.
- Added original-name compatibility wrappers.
- Added debug and optimized strict validation workflows, four numerical test
  suites, a demonstration, an integration example, and a dated CSV fitter.
- Preserved GPL-2.0-or-later licensing in `LICENSE`, `fpm.toml`, and every
  Fortran source file.
