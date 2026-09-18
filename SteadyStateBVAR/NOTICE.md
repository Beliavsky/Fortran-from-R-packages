# Notices and attribution

This directory is an unofficial modern Fortran translation of computational code from **SteadyStateBVAR 0.2.0**.

Upstream package:

- Name: SteadyStateBVAR
- Version: 0.2.0
- Author, maintainer, and copyright holder listed by the upstream package: Mark Becker
- Upstream source: https://github.com/markjwbecker/SteadyStateBVAR
- Upstream package license: GNU General Public License, version 3 or later
- Upstream citation: Becker, M. (2026), *SteadyStateBVAR: Bayesian Vector Autoregressions with Steady-State Priors*.

The upstream GPL notice from `inst/stan/include/license.stan`, the upstream `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, `CITATION`, the computational R sources used for this translation, and the four Stan model sources are retained under `upstream/` for provenance. The Fortran translation is distributed under GPL-3.0-or-later; see `LICENSE`.

This translation reuses existing top-level projects from `Beliavsky/Fortran-from-R-packages` through FPM sibling path dependencies rather than copying their sources:

- `rfortran-core` (MIT): `r_kinds` for `dp`, `r_distributions` for `r_qnorm`, and `r_quantiles` for R type-7 quantiles.
- `rfortran-linalg` (MIT): checked square/SPD solves, Cholesky factors, matrix inversion, and SPD inverse/log-determinant operations.
- `MTS` / FPM package `mts-fortran` (Artistic-2.0): `var_psi_weights`, the translated counterpart of `MTS::VARpsi` used by conditional forecasts and IRFs.

No dependency source is vendored in this package. `rfortran-linalg` owns the provenance and license obligations of its LAPACK backend. This package does not call LAPACK directly and therefore does not declare a direct `fortran-lapack` dependency. Federico Perini is the original author of the `fortran-lapack` project used by the shared linear-algebra layer; see the `rfortran-linalg` package for its pinned backend revision and full attribution.

The translation is not endorsed by Mark Becker, CRAN, the R Foundation, Stan, MTS, or the authors of the shared Fortran dependencies.
