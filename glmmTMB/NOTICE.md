# NOTICE and provenance

## Upstream project

- Package: `glmmTMB`
- Upstream version: `1.1.14`
- Upstream project: `https://github.com/glmmTMB/glmmTMB`
- Upstream license declaration: `AGPL-3`
- Supplied archive used for this translation: `glmmTMB-master.zip`
- Upstream package date in `DESCRIPTION`: 2026-01-14

The original author and contributor list is preserved verbatim in
`upstream/DESCRIPTION`.  The upstream citation metadata is preserved verbatim in
`upstream/CITATION`; it includes the Brooks et al. (2017) R Journal article and
the McGillycuddy et al. (2025) Journal of Statistical Software article.

## Computational sources translated

The portable Fortran implementation was derived principally from these files in
the supplied upstream archive:

- `src/glmmTMB.cpp` — family/link codes, conditional likelihood, zero inflation,
  zero truncation, Gaussian random-effect covariance structures, and priors;
- `src/distrib.h` — beta-binomial, generalized Poisson, skew-normal, Cauchy,
  Bell/Lambert-W, COM-Poisson variance support, and distribution helpers;
- `src/cordistrib.h` — LKJ, Wishart, inverse-Wishart, and multivariate-gamma
  correlation/covariance prior helpers;
- `R/family.R` — family variance functions and family-level computational
  parameterizations.

The supplied upstream `MD5` file is copied to `upstream/MD5` so the exact source
snapshot can be cross-checked.

## Reused dependency

The translation uses the sibling FPM path dependency `../TMB`, which is the
modern Fortran translation of TMB prepared for the same repository.  It supplies
explicit-interface normal/beta/gamma/binomial/Student-t distribution kernels,
multivariate-normal likelihood evaluation, Cholesky factorization, and the
unstructured-correlation transformation.

No TMB, BLAS, LAPACK, Eigen, RcppEigen, or translated R-package source is copied
or vendored into this directory.

## Translation changes

This is a source-language translation, not a verbatim source copy.  C++
templates and TMB AD types are represented by `real(dp)` numerical APIs.  R and
TMB simulation functions that obtain random numbers from R are excluded.
The arbitrary-order Matérn Bessel K calculation is implemented with numerical
quadrature.  The Tweedie density is evaluated from its compound-Poisson/gamma
series rather than by importing another R-package translation.

No warranty is provided.  See `LICENSE` for the governing terms.
