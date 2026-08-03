# Porting notes

## Data representation

R matrices and data frames are replaced by rank-2 `real(dp)` arrays. Rows are
observations and columns are variables. Logical input is not overloaded; callers
can convert logical values to zero/one real values explicitly.

## Association matrices

Pearson calculations use an ordinary sample correlation matrix. Spearman first
replaces each column with average ranks and then uses Pearson correlation.
Kendall computes pairwise tau-b, including tie corrections.

The upstream Pearson path uses a covariance matrix. Both partial and
semi-partial correlation formulas are invariant to diagonal positive rescaling,
so covariance and correlation produce the same coefficients.

## Generalized inverse

The upstream package calls `MASS::ginv()` when the covariance determinant is
near zero. This port symmetrizes the association matrix, computes a Jacobi
eigendecomposition, and inverts eigenvalues above a relative threshold. This is
more stable than using the determinant as a singularity test. The retained rank
and fallback flag are returned to the caller.

## Significance tests

Pearson and Spearman use the same Student-t statistic and residual degrees of
freedom as ppcor 1.1. Kendall uses the upstream large-sample normal statistic.
The Student-t CDF is evaluated through the regularized incomplete beta function.

## Differences from R behavior

- Non-finite values are rejected with `ppcor_nonfinite_data`; there is no `NA`
  propagation or pairwise-complete mode.
- A constant column returns `ppcor_constant_variable`.
- Degenerate coefficients are clipped to `[-1,1]` before significance testing,
  avoiding NaNs from tiny pseudoinverse roundoff.
- R warnings are represented by `result%message` and
  `result%used_pseudoinverse`.
- R vector recycling, names, dimnames, and print methods are omitted.
