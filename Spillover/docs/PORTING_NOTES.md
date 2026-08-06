# Porting notes

## Numerical conventions

The fitted VAR is written as

```text
y(t) = c + d*t + A(1)y(t-1) + ... + A(p)y(t-p) + e(t).
```

All equations use the same regressor matrix and are estimated by ordinary least squares. The residual covariance divides the residual cross-product by `n_effective - ncoef`, matching the degrees-of-freedom adjustment in `summary.varest`.

The generalized FEVD uses

```text
theta(i,j;H) = sigma(j,j)^(-1) * sum_h (e_i' Phi_h Sigma e_j)^2
               / sum_h e_i' Phi_h Sigma Phi_h' e_i.
```

When normalization is requested, each row is divided by its row sum. Orthogonalized FEVDs use the lower Cholesky factor of the innovation covariance.

## Ordering averages

The R package delegates ordering averages to `fastSOM`. This port is independent:

- sampled averages use a deterministic Park-Miller generator and Fisher-Yates permutations;
- exact averages enumerate permutations lexicographically;
- exact enumeration defaults to dimensions no larger than nine and can be changed with `exact_limit`.

Consequently sampled results are reproducible within this library but do not share R's random-number stream. Exact values are deterministic.

## Upstream inconsistencies

Two source-compatible options are provided:

1. In `O.spillover.R`, `partial` calls the exact routine and `total` calls the estimated routine, contrary to both the documentation and `roll.net.R`. The Fortran default follows the documentation. Set `source_compatible=.true.` in `orthogonalized_spillover` or `o_spillover` to reproduce the source labels.
2. In `dynamic.spillover.R`, the quantities named `from` and `to` are formed from column and row sums respectively, opposite the package's own table definitions. The Fortran default uses row sums for `from` and column sums for `to`. Set `source_compatible_direction=.true.` to reproduce the R assignment.

## Deliberate omissions

- dates and `zoo` indexes;
- missing-value filling and R's `na.omit` behavior;
- formula parsing, seasonal dummies, and exogenous regressors from `vars::VAR`;
- plotting and long-form data-frame construction;
- parallel execution.

Inputs must be finite numeric arrays. Failures are returned through integer status codes and optional diagnostic messages.
