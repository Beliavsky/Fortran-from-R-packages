# Provenance and license notices

This package translates the computational routines of **trend 1.1.8**, authored
and maintained by **Thorsten Pohlert**. Upstream release date: 2026-09-03.

- CRAN: https://cran.r-project.org/package=trend
- Source mirror: https://github.com/cran/trend
- Tested source revision: `ee55a676a23b2d080d953fe82e627216562f45bf`
- Target repository inspected at revision:
  `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`

The preserved upstream headers state Copyright (C) 2015–2018 Thorsten Pohlert
(with individual files specifying 2016–2018, 2017–2018, or 2018). These headers
are retained verbatim in `upstream/`. Upstream DESCRIPTION declares GPL-3;
the computational source headers expressly permit GNU GPL version 3 **or any
later version**. The new `itam.R` file has no broader per-file license grant. The combined
translation conservatively uses GPL-3.0-only, retaining all upstream notices
that permit later versions for individual source files.
The complete GNU GPL version 3 is in `LICENSE`.

This software is provided WITHOUT ANY WARRANTY, including implied warranties
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE, as detailed in LICENSE.

## Preserved reference source

`upstream/DESCRIPTION`, `upstream/NAMESPACE`, computational R source files,
internal `utilfn.R`, and original C/fixed-form Fortran Monte Carlo source files
are unchanged reference material. FPM does not compile this directory. The
maintained translation in `src/`, tests, and examples is free-form Fortran.
Plotting, presentation methods, and package datasets are not translated. The
original plot.itam function is retained only as part of the unchanged itam.R
reference file.

## Changes in this translation

The R/S3/time-series interface becomes explicit arrays, periods, result records,
and integer diagnostic statuses. Algorithms are rewritten as pure, explicit-
interface free-form Fortran using the shared `r_kinds:dp`. Numerical primitives
and compatible translated algorithms are linked through sibling FPM dependencies.
Monte Carlo uses a local Park-Miller/Box-Muller generator or caller-provided null
samples rather than R's process-wide RNG. R's observed sample standard deviation
and simulated population standard deviation conventions are preserved. Numerical
validation and documented edge checks have been added. Centering/scaling corrects
raw-moment cancellation in WW and stabilizes change-point statistics. No claim of endorsement
by the upstream author is made.

## Dependency reuse and attribution

No dependency source is copied into this package.

- `../rfortran-core` (MIT): `dp`, ranking, sorting, median, mean, sample standard
  deviation, correlation, and normal/Student-t probabilities and normal quantiles.
- `../Kendall` (GPL-2.0-or-later), A.I. McLeod: compatible `mann_kendall`
  score, tie-corrected variance, and tau-b. The trend adapter supplies its own
  selected-tail significance and continuity option. The two-observation seasonal
  special case is handled directly because the shared interface requires three.
- `../randtests` (GPL-2.0-or-later), Frederico Caeiro and Ayana Mateus:
  compatible `bartels_rank_test` RVN and lower-tail beta probability. Trend's
  distinct normal approximation for n >= 100 is applied by the adapter.
- `../exactRankTests` (GPL-2.0-or-later), Torsten Hothorn and Kurt Hornik:
  `wilcox_exact` rank-sum computation and exact/normal significance for Lanzante.
  The adapter restores `stats::wilcox.test`'s tie-based default exact selection
  and normal continuity correction.

The upstream `extraDistr::pnsbeta(T,a,a,min=0,max=4)` call is the beta CDF at
`T/4`. The existing randtests implementation performs precisely this transform
using core regularized beta, and its result is checked against actual trend R.
The existing extraDistr translation is GPL-2.0-only, so it is not linked into this
GPL-3 translation. Existing randtests Cox-Stuart and stats rank-test interfaces
were inspected: randtests uses half-series binomial comparisons (trend uses
thirds and a normal approximation), while stats lacks the required exact
Wilcoxon option. These are not interchangeable.

The Park-Miller recurrence and Box-Muller equations are mathematical algorithms,
implemented here without copying a third-party RNG implementation. References:
Park and Miller (1988), Communications of the ACM 31(10), 1192–1201;
Box and Muller (1958), Annals of Mathematical Statistics 29(2), 610–611.
