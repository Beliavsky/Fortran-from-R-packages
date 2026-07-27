# Porting notes

## Data orientation

As in HDShOP, a return matrix is `p` by `n`: assets are rows and observations are
columns.

## Linear algebra

The translation includes a pivoted Gauss-Jordan inverse and a symmetric Jacobi
eigensolver. High-dimensional `p > n` portfolio routines use an eigenvalue-based
Moore-Penrose inverse, replacing `MASS::ginv`. Nonlinear Ledoit-Wolf shrinkage
uses the same analytical kernel formulas and ascending eigenvalue ordering as
the R implementation.

## Missing observations

The original sample-covariance routine computes row means with `na.rm=TRUE` but
then performs an unguarded cross-product, which can propagate missing values.
The Fortran routine uses pairwise finite observations for each covariance entry.
For fully observed arrays, its result is identical to the R formula.

## Typed results

R lists and S3 classes are represented by `matrix_shrink_result`,
`mean_shrink_result`, `portfolio_result`, `mvsp_test_result`, and
`frontier_result`. Each result has an `ok` flag and diagnostic message.

## High-dimensional behavior

The `p < n` routines return asymptotic per-weight intervals and test statistics.
The upstream package omits these intervals for `p > n`; the Fortran translation
does the same. Square `p = n` inputs are rejected by the shrinkage dispatcher,
matching the upstream behavior.

## Numerical safeguards

- Explicit covariance symmetrization avoids compiler-dependent asymmetry.
- Near-zero eigenvalues are excluded when forming pseudoinverses.
- Negative roundoff in variances is clipped before square roots.
- Normal and chi-square-one probabilities are implemented internally.
- Infinite risk aversion is handled as an exact zero reciprocal, producing GMV
  weights without arithmetic involving infinity.

## Random covariance matrices

`random_covariance_matrix` follows the upstream Wishart-eigenvector construction
and supports deterministic seeds. The resulting matrix has the requested
spectrum up to eigensolver tolerance.
