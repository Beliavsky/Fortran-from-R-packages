# Porting notes

## Design

The R package is object-oriented and delegates several estimators to
`robustbase` and `pcaPP`. The Fortran port uses plain derived types and explicit
subroutines. Data matrices use observation-major mathematical indexing
`x(n, p)`; Fortran still stores these arrays column-major internally.

## Robust covariance details

- MCD uses random elemental starts, repeated h-subset concentration steps,
  median chi-square consistency scaling, and 97.5% chi-square reweighting.
- MVE evaluates elemental ellipsoids and h-subsets using a log-volume objective.
- OGK follows the orthogonalization/reconstruction structure of `covOPW.c` and
  supports MAD or tau scales and GK or quadrant pairwise association.
- S and MM estimators use Tukey-biweight redescending weights with robust MCD
  initialization. Upstream tuning-constant integration code is not reproduced
  exactly; stable dimension-aware chi-square tuning is used.
- SDE constructs coordinate, pair-difference, and random unit directions, then
  forms a high-breakdown subset from maximum standardized projection distances.
- MRCD uses a robust diagonal Qn target and increases regularization until the
  covariance condition number is controlled, unless the caller supplies `rho`.

## PCA details

`pca_hubert` performs classical dimension reduction followed by MCD in score
space and rotates the loading basis by the robust score covariance. This covers
the core robust subspace logic without reproducing R plotting, skew-adjusted
outlier maps, or every pcaPP initialization branch.

`pca_grid` and `pca_proj` search coordinate and random orthogonal directions and
maximize Qn scale. This is a self-contained projection-pursuit implementation.

## LdaPP and Linda

The R implementations contain extensive formula interfaces and method-specific
robust univariate optimization. The Fortran `lda_pp_fit` exposes a robust
projection/covariance discriminant fit, while `linda_fit` uses MCD group centers
and a robust common residual covariance. Prediction formulas match standard
LDA.

## Error handling

R exceptions are represented by integer status fields. The public constants are
`rrcov_success`, `rrcov_invalid_argument`, `rrcov_dimension_error`,
`rrcov_singular`, `rrcov_no_convergence`, and `rrcov_allocation_error`.

## Reproducibility

Randomized estimators accept an integer `seed`. The internal generator is the
Fortran processor generator seeded deterministically; it does not reproduce R's
Mersenne-Twister stream.
