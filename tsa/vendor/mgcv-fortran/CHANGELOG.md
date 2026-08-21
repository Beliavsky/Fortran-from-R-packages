# Changelog

## 0.1.1 - 2026-08-05

- Fixed the bundled FPM dependency manifest: the dependency key `splines`
  now matches `dependencies/splines/fpm.toml` package name `splines`.
- This resolves FPM's `Dependency name 'splines-fortran' found, but expected
  'splines' instead` model error.
- No Fortran source or numerical behavior changed.

## 0.1.0 - 2026-08-05

- Added a modern Fortran 2018/FPM numerical core for dense generalized
  additive models.
- Added penalized IRLS and smoothing-parameter selection by GCV, UBRE, and a
  REML-like criterion.
- Added cubic regression, P-spline, cyclic, low-rank radial, random-effect,
  and tensor-product smooth construction.
- Added dense matrix, distribution, simulation, constrained-fit, and
  discrete cross-product utilities.
- Bundled the translated `splines` package as a local FPM dependency.
- Preserved GPL-2.0-or-later licensing and original source provenance.
- Omitted plotting and R-specific formula/S3/data-frame infrastructure.
