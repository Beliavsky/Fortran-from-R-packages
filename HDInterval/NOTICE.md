# Notices and attribution

This package is a numerical Fortran translation of **HDInterval 0.2.4**,
“Highest (Posterior) Density Intervals”, dated 2022-11-16, published on CRAN
2022-11-17. The upstream `DESCRIPTION` identifies:

- Mike Meredith, author
- John Kruschke, author
- Ngumbang Juat, creator/maintainer (ngumbangjuat@hotmail.com)

Upstream source is licensed **GPL-3**. This translation and its newly maintained
Fortran, tests, example and accompanying original documentation are distributed
under **GPL-3.0-only**. A full GNU General Public License version 3 is in
`LICENSE`; upstream's GPL-3 designation has not been broadened to “or later”.
Unmodified upstream files retain their original notices and attribution.

The upstream function-method comments credit Kruschke (2011), *Doing Bayesian
Data Analysis*, for the inverse-CDF minimum-width interval construction, and
mention the related TeachingDemos package by Greg Snow. Those comments are
preserved in `upstream/R/hdi_function.R`. This is not a claim that either author
reviewed or endorsed the Fortran translation.

## Shared dependency and compatibility review

`rfortran-core` is used as an unmodified sibling FPM path dependency for `dp` and
NaN-aware sorting. Its repository license is MIT; source is not bundled here.

The sibling `stats` numerical module was inspected before implementing scalar
search logic. It offers pure `optimize` and `uniroot` operations with apparently
suitable APIs, but an independent FPM smoke build of the unchanged package failed
because `r_stats_multivariate` imports `r_matrix_statistics`, for which this
repository checkout has no module definition. Its transitive numerical dependency
was resolved through the Beliavsky fortran-lapack fork in that smoke test, so the
failure was not a missing system BLAS/LAPACK library. `stats` could not be reused
as a functioning FPM dependency without modifying a sibling package.

Accordingly, only focused HDInterval-specific golden-section interval-width
search and expanding CDF bisection are implemented in `hdinterval.f90`. No
`stats`, `r.f90`, `r_mod.f90`, BLAS, LAPACK, ARPACK or translated dependency source
has been copied or embedded. No shared sibling source was changed. No linear
algebra dependency or external system-library link is required by this package.

Fortran translation and validation prepared in 2026. See `PROVENANCE.md` for
source integrity information and `API_COVERAGE.md` for the exact scope.
