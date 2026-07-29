# Porting notes

## Model convention

For an ordinary normal variance mixture,

```text
X = location + sqrt(W) A Z,
```

where `A A^T = scale` and `Z` is standard multivariate normal.

For grouped mixtures, coordinate `j` uses `W(groupings(j))`. Conditional on the
mixing vector, the covariance is therefore

```text
D(sqrt(W_group)) scale D(sqrt(W_group)).
```

The Fortran `nvmix_model` stores location, scale, group assignments, mixing
families and mixing parameters explicitly.

## Mixing families

- `mix_constant`: `W = 1`.
- `mix_inverse_gamma`: `W ~ InvGamma(df/2, df/2)`.
- `mix_pareto`: quantile `(1-u)^(-1/alpha)`.
- `mix_gamma`: mean-one gamma with shape `a` and scale `1/a`.

## Integration

Closed-form normal and ungrouped Student densities are used whenever possible.
Other densities integrate the conditional normal density over deterministic
Halton points and use log-sum-exp stabilization.

Rectangle probabilities are estimated from deterministic quasi-random mixture
and normal points. Independent batches provide a standard-error estimate.
This is simpler than the upstream preconditioned adaptive RQMC/C implementation,
but keeps the model definition and target integral unchanged.

## Student fitting

For each candidate degrees of freedom, the location and scale are updated with
the standard multivariate-t EM weights

```text
w_i = (df + d) / (df + Mahalanobis_i^2).
```

A bounded golden-section profile search estimates `df`. Fixed-`df` fitting is
also supported.

## Copula fitting

The Student-copula estimator profiles over degrees of freedom, transforms each
margin with the t quantile, and estimates the correlation matrix from the
transformed observations. The grouped estimator estimates marginal tail
parameters and then the transformed correlation. These are intentionally
self-contained and do not duplicate all upstream ECME or third-party optimizer
options.

## Skew-t

The translated skew-t follows the variance-mean-mixture representation

```text
X = location + gamma W + sqrt(W) A Z,
W ~ InvGamma(df/2, df/2).
```

Density and distribution calculations integrate over `W`; simulation uses the
same representation directly.

## Dependence measures

For grouped mixtures, covariance is computed from `E(W)` and `E(sqrt(W))`.
The standard elliptical Kendall formula is exact for ordinary normal variance
mixtures. The equal-df Student tail-dependence formula is exact. Unequal grouped
degrees of freedom use a documented harmonic-mean approximation rather than the
upstream RQMC integral.

## Numerical safeguards

- Positive-definite scale matrices are checked by Cholesky factorization.
- Density integration uses log-sum-exp stabilization.
- Probabilities are clipped to `[0,1]` after numerical integration.
- Quantile solvers expand brackets before bisection.
- Random seeds are explicit and deterministic.
- All source avoids assumptions about short-circuit logical evaluation.
