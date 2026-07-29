# Porting notes

## Data orientation

The Fortran predictor matrix has shape `(n_observations, n_predictors)`. The
intercept is inserted internally and cannot be removed.

## Closed-form transformation

The upstream estimator fits

```text
Y_t = beta_0 + beta' u_t + epsilon_t,
Y_t = Phi^-1(loss_t).
```

With `sigma2 = mean(epsilon_t^2)`, the original-space parameters are

```text
p       = Phi(beta_0 / sqrt(1 + sigma2))
rho     = sigma2 / (1 + sigma2)
kappa_j = beta_j / sqrt(1 + sigma2).
```

These equations are preserved directly.

## Finite-portfolio correction

When `portfolio_size` is present, observed rates in `[0,1]` are shrunk around
their sample mean using the Yang correction before the probit transformation.
The implementation rejects zero sample variance and non-positive corrected
variance explicitly.

## Covariance estimation

The IID covariance reproduces the upstream construction:

- conventional OLS covariance for regression coefficients;
- the upstream variance formula for `sigma2`;
- zero cross-covariance between OLS coefficients and `sigma2`;
- the same analytical delta-method Jacobian for `(p, rho, kappa)`.

The HAC path reproduces the upstream influence functions but uses a native
Bartlett/Newey-West long-run covariance instead of `sandwich::lrvar`. The
result is divided by the observation count to estimate the covariance of the
sample estimator.

## Numerical choices

- The normal inverse CDF uses a rational approximation followed by Newton
  refinement.
- OLS uses normal equations with pivoted matrix inversion. This is appropriate
  for the package's small macro-factor designs; nearly singular designs are
  rejected.
- Covariance matrices are explicitly symmetrized so Windows and Linux builds
  return identical opposite-triangle entries.
- All public computations use `real(dp)`, where `dp = kind(1.0d0)`.

## Random generation

`random_vasicek` uses Fortran's intrinsic pseudorandom generator after filling
its seed vector deterministically from an optional 64-bit integer. A fixed seed
is reproducible within a compiler/runtime, but the exact stream is not promised
to match R or other Fortran runtimes.
