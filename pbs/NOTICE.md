# NOTICE and provenance

## Upstream package

- Package: `pbs`
- Version: 1.1
- Title: *Periodic B Splines*
- Author/Maintainer: Shuangcai Wang
- Upstream date: 2013-03-22
- License: GPL-2
- Upstream repository metadata: CRAN

The original package metadata, namespace, R source, and manual page are preserved under `upstream/`.

The upstream implementation states that `pbs()` is modified from the B-spline functionality in R's `splines` package and uses `spline.des()`/`splineDesign()` for spline evaluation. This Fortran translation independently implements the corresponding normalized B-spline basis and derivative recurrence rather than copying R's implementation.

## Translation provenance

This translation was prepared for the `Beliavsky/Fortran-from-R-packages` repository. It is a source translation/independent reimplementation of the package's numerical behavior and remains subject to the upstream GPL-2 licensing terms. See `LICENSE`.

No BLAS, LAPACK, R runtime source, `r.f90`, `r_mod.f90`, or translated R-package dependency is copied into this package.
