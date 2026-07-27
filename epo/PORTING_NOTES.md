# Porting notes

## Array orientation

R stores a return data frame or matrix with observations in rows and assets in
columns. The Fortran API uses the same logical orientation:

```text
returns(number_of_observations, number_of_assets)
```

Fortran is column-major internally, but this does not alter the interface.

## Correlation shrinkage

The upstream calculation is

```text
shrunk_cor = (1 - w) correlation + w identity
cov_tilde = standard_deviation_matrix * shrunk_cor *
            standard_deviation_matrix
```

The implementation also uses the equivalent identity

```text
cov_tilde = (1 - w) covariance + w diagonal(covariance).
```

Both `shrunk_correlation` and `shrunk_covariance` are returned.

## Linear algebra

The R code calls `solve(cov_tilde)`. The Fortran code avoids forming an
explicit inverse and solves the positive-definite system using Cholesky
factorization. A singular or non-positive-definite matrix returns a typed
failure rather than unstable weights.

## Endogenous anchored scaling

The upstream denominator is written as

```text
signal' inverse(C_tilde) C_tilde inverse(C_tilde) signal.
```

For symmetric nonsingular `C_tilde`, this equals

```text
signal' inverse(C_tilde) signal.
```

The Fortran implementation uses the simplified expression after solving the
linear system. This changes neither the mathematics nor valid numerical
results and avoids two explicit matrix inversions.

## Validation differences

The R package documents `w` as lying between zero and one but checks only that
it is numeric. The Fortran port enforces the documented interval.

The R package allows any numeric `lambda`; the Fortran port requires a positive
value whenever `lambda` is used. In endogenous anchored EPO, `lambda` is not
used and therefore is not required to be positive.

R's default covariance calculation propagates missing values. The Fortran port
rejects non-finite values explicitly rather than returning a partially invalid
portfolio.

## Normalization

Normalization divides all weights by their sum. It is not a constrained
optimization and does not enforce nonnegative holdings. If the raw weights sum
to zero within a scale-aware floating-point tolerance, the routine reports a
normalization failure.

## Full shrinkage

For simple EPO, `w = 1` removes all correlations but still uses the asset
variances and signals. It does not generally produce equal weights.

For anchored EPO, `w = 1` returns the anchor exactly before normalization,
matching the upstream implementation.
