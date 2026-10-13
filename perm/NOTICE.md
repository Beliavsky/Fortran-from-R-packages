# Notices and provenance

## Upstream

- Package: **perm**, version **1.0-0.4**, dated 2023-08-24.
- Author/maintainer: **Michael P. Fay**, mfay@niaid.nih.gov.
- Official source: <https://cran.r-project.org/src/contrib/perm_1.0-0.4.tar.gz>.
- Official package page: <https://cran.r-project.org/package=perm>.
- Source tarball SHA-256: `561e01793a9032b6a3676a7819fe24175eeb9e96a5314bd78b9c48d46e627be9`.
- Upstream DESCRIPTION records `License: GPL`; CRAN identifies the permitted
  choice as GPL-2 or GPL-3. Both complete texts are included, and the translation
  retains that choice (`GPL-2.0-only OR GPL-3.0-only`).
- `upstream/DESCRIPTION`, `upstream/NAMESPACE`, `upstream/MD5`, original R files,
  and manual pages are retained unchanged as provenance. The MD5 file describes
  the original archive; files omitted from this translation are not asserted
  to be present. No upstream third-party dependency code is embedded.
- The package's reference is Fay and Shaw, “Exact and Asymptotic Weighted
  Logrank Tests for Interval Censored Data: The interval R Package,” Journal of
  Statistical Software 36(2), <https://doi.org/10.18637/jss.v036.i02>.
- The Monte Carlo add-one correction follows the upstream citation to Fay,
  Kim and Hachey (2007), JCGS, pages 946–967, equation 5.3.
- The network algorithm's upstream reference is Agresti, Mehta and Patel
  (1990), JASA, pages 453–458.

## Translation

Maintained Fortran, test programs, examples and new documentation are a 2026
translation by the repository's Fortran translation contributors. Numerical
operations are adapted from the R files named in `fpm.toml`. Network paths are
represented by grouped-score occupancy counts rather than R lists and graph
edge vectors; each path has the same product of binomial multiplicities.
Log-weight normalization avoids overflow of `choose(n,m)`. No R formula,
printing, warning, S3-object or dynamic callback implementation is claimed.

## Sibling reuse and compatibility decisions

Repository inspected at commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24` before implementation.

- `rfortran-core` is a sibling FPM path dependency, not a copied component. Its
  MIT license carries Copyright (c) 2026 Beliavsky. The translation imports
  `dp` from `r_kinds`, sorting from `r_ordering`, distribution probabilities
  and the F quantile from `r_distributions`, regularized beta and log-choose
  from `r_special`, and stable `expm1`/`log1p` from `r_transforms`.
- `exactRankTests` was inspected: its distribution API maps arbitrary scores
  onto a rounded lattice. That is incompatible with perm's exact arbitrary
  real-score occupancy enumeration, so it is not used or copied.
- `extraDistr` was inspected for beta quantiles. Its inverse is impure and uses
  an absolute probability-scale tolerance; the focused core provides the
  required CDF/F-quantile chain. The pCI calculation uses a beta/F transform,
  analytic shape-one boundaries and relative-tolerance refinement against the
  shared regularized beta CDF. No extraDistr source is reproduced.
- No `stats`, compatibility umbrella, BLAS, LAPACK or ARPACK dependency is
  needed. There are no system-library links or external numerical binaries.

## Reference verification

`reference/generate_reference.R` was executed with the official installed
perm 1.0-0.4 under R 4.5.0. `reference/r_oracle.txt` records the resulting
numerical fixtures. It is generated data, not a substitute for the retained
source. The Fortran test additionally evaluates independent tiny enumerations
and deterministic Monte Carlo calibration cases.
