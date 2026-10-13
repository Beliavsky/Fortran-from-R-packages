# License and upstream provenance

This package is a modern free-form Fortran translation of the computational code of
R package **PropCIs 0.3-0**, by **Ralph Scherer**, maintainer
`shearer.ra76@gmail.com`. Upstream date: 2018-02-22; CRAN publication:
2018-02-23 16:49:49 UTC.

- Official source: https://cran.r-project.org/src/contrib/PropCIs_0.3-0.tar.gz
- Upstream repository: https://github.com/shearer/PropCIs
- Official source archive SHA-256:
  `cd35775f4d36e642663e727450c53708f17b3c4340e1bc2c1752fd17118a9ffb`
- Verified source retrieved 2026-10-11.
- `upstream/DESCRIPTION`, `upstream/NAMESPACE`, and `upstream/R/*.R` are unmodified
  copies of this package's source and metadata for attribution and audit.
- Maintained Fortran code, tests, examples, documentation, and reference tooling
  were adapted in 2026 by Fortran translation contributors using OpenAI.

Upstream declares `License: GPL`, without specifying a version. This distribution
elects **GPL version 2**, as permitted by the unversioned grant and GPL section 9.
The full license is in LICENSE. Retained upstream source and metadata have not had
their original license declaration rewritten. No additional upstream copyright or
license files were present in the official package archive.

The translation preserves the interval definitions and attribution. Numerical
methods are intentionally changed where documented in README.md and API_COVERAGE.md:
bounded constrained-score roots, stable odds quadratics, deterministic posterior
quadrature/inversion, connected-interval width optimization, explicit statuses,
and fixes for verified upstream defects. These changes are not upstream releases.

## Reused repository components

All dependencies are separate sibling FPM packages; none of their implementation
sources, notices, licenses, BLAS, LAPACK, or R compatibility copies are embedded here.
Their notices/licenses remain with those dependencies:

- `rfortran-core`, Vivek Rao and contributors, MIT: kind, normal quantile, beta,
  digamma, trigamma, and other special functions.
- `extraDistr`, upstream extraDistr authors and Fortran translation contributors,
  GPL-2.0-only: existing beta quantile and binomial helpers.
- `RcppNumerical`, Yixuan Qiu and contributors, with the repository's OpenAI/Fortran
  translation, GPL-2.0-or-later: adaptive Gauss-Kronrod integration.
- `lbfgs`, upstream lbfgs/libLBFGS copyright holders and translation contributors,
  GPL-2.0-or-later; `lbfgsb3`, upstream lbfgsb3c authors and translation contributors,
  GPL-2.0-only: declared transitive dependencies of RcppNumerical.

The selected GPL-2 distribution is compatible with these dependency licenses.
No system BLAS or LAPACK libraries are linked. This package does not require
fortran-lapack, so no fork or LAPACK implementation is incorporated.
