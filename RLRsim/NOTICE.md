# NOTICE and provenance

This is a clean Fortran translation of computational code from RLRsim version 3.1-9.

Upstream project: RLRsim
Upstream URL: https://github.com/fabian-s/RLRsim
Upstream CRAN package version: 3.1-9
Upstream primary author and maintainer: Fabian Scheipl
Upstream contributor listed in DESCRIPTION: Ben Bolker
Upstream license declaration: GPL

The original package metadata, NAMESPACE, NEWS, CITATION, R sources relevant to coverage, and the original C++ computational kernel are retained under `upstream/` for attribution and traceability. No upstream Rcpp, R, BLAS, LAPACK, or translated dependency source is embedded in `src/`.

The Fortran translation preserves the numerical structure of the R/C++ simulation algorithm while replacing R-specific object handling and RNG facilities with explicit Fortran APIs. See `README.md` and `API_COVERAGE.md` for compatibility differences.

The upstream DESCRIPTION declares `License: GPL`. Copies of GPL version 2 and GPL version 3 are included in `LICENSES/` to preserve the applicable GPL license texts used by R package distributions. This translation is distributed under GPL-2.0-or-later in `fpm.toml`.

No direct LAPACK dependency is used, so the Federico Perini / Beliavsky `fortran-lapack` fork is not required by this package.
