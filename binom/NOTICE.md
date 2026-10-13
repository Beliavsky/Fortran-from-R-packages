# Attribution and provenance

This package translates the computational code of CRAN **binom 1.1-2**.
Upstream authors: **Sundar Dorai-Raj** and **Arianna D. Cascone** (current
maintainer). Original package title: “Binomial Confidence Intervals for Several
Parameterizations”. Preserve their authorship and all original source notices.

- Official archive: https://cran.r-project.org/src/contrib/binom_1.1-2.tar.gz
- CRAN publication: 2026-07-07; DESCRIPTION date: 2026-07-06.
- SHA-256: `375f7c080a71bc259b0e1fc69b5cfbefc9a4724881d2fdd29681236952c878ea`.
- Original archive MD5 manifest: all 35 listed files verified on extraction.
- Original R and C computational sources, DESCRIPTION, NAMESPACE, MD5 and NEWS
  are retained under `upstream/` as provenance only. That directory is outside
  FPM's library source tree. No R or C runtime is required.

Upstream DESCRIPTION declares **GPL** without a version restriction. This
translation is offered under GPL-2.0-or-later. The full GPL v2 text is `LICENSE`;
GPL v3 is also included in `LICENSES/GPL-3.0.txt`. The original license declaration
is retained unchanged in `upstream/DESCRIPTION`. No separate upstream copyright
notice was present in the translated R and C source files.

## Dependencies and licensing of combined programs

Dependencies are referenced, not bundled. Repository inspection used
Beliavsky/Fortran-from-R-packages commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.

- `rfortran-core`: Vivek Rao and contributors, MIT. Shared real kind, normal
  distributions, beta CDF/log beta, stable sorting and elementary transformations.
- `rfortran-linalg`: Vivek Rao and contributors, MIT. Pure high-level
  `solve_system` solves the small FMM profile interpolation system. Its existing
  tested transitive fortran-lapack pin is retained unchanged.
- `fortran-lapack`: Federico Perini and contributors. Used transitively through
  rfortran-linalg; its own notices and licenses apply. No direct LAPACK calls or
  separately installed BLAS/LAPACK are required. The Beliavsky fork is the
  repository's preferred route for direct LAPACK/pure call chains; this package
  instead uses the compatible existing shared high-level API and does not
  override its older tested dependency.
- `lbfgsb3`: existing top-level Fortran translation of lbfgsb3c. Its wrapper is
  GPL-2.0-only and the numerical kernel BSD-3-Clause. See that sibling package's
  NOTICE.md and LICENSES for Matthew L. Fidler, John C. Nash, Ciyou Zhu, Richard Byrd, Jorge Nocedal,
  Jose Luis Morales and collaborators. A linked program including this optimizer
  must comply with GPL version 2; the translation's “or later” option does not
  relicense the dependency.

## Focused numerical code and compatibility inspection

The existing `splines` sibling offers natural cubic interpolation, not the FMM
endpoint conditions used by R's default `stats::spline`. The existing `stats`
translation has an FMM operation but its full FPM package requires an unavailable
`r_matrix_statistics` module, and its interpolation API is not an independently
usable dependency. Neither package is copied or vendored. `binom_profile.f90`
implements only the small profile-specific FMM assembly/evaluation, using the
shared linear solver. Its endpoint third derivatives are the cubic divided
differences of the first/last four knots (Forsythe, Malcolm and Moler, 1977),
matching R's `stats::spline(method="fmm")` convention. This is not a copy of the
existing translated stats dependency module.

The core library provides the regularized beta CDF but no beta quantile. The
local tail-preserving beta inversion is therefore a small numerical adapter,
not a copied translated package. Highest-density beta intervals solve the same
posterior-mass/equal-density equations as upstream's C implementation using a
bounded allocation of the excluded tail probability. The Brent root iteration is translated from the retained upstream binom.c,
with the actual R uniroot default tolerance preserved for LRT. All inputs and numerical
changes are documented in README.md and API-COVERAGE.md.
