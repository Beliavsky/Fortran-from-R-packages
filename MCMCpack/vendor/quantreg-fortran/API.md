# API

All real data use `real(dp)` with `dp = kind(1.0d0)`.

## Result types

- `rq_result`: coefficients, residuals, dual values, tau, iteration counters,
  and an integer `info` status.
- `rq_multi_result`: coefficient matrix for multiple quantiles.
- `lprq_result`: evaluation grid, fitted local quantile and local derivative.
- `nlrq_result`: nonlinear coefficients, residuals and objective.

## Dense quantile regression

### `rq_fit_fnb(x, y, tau, result [, beta, eps])`
Frisch-Newton interior-point QR. `x` is `(n,p)`.

### `rq_fit_fnc(x, y, R, r, tau, result [, beta, eps])`
Inequality-constrained QR with the same convention as upstream:
`matmul(R,b) >= r`.

### `rq_fit_qfnb(x, y, taus, result)`
Fits a vector of quantiles.

### `rq_wfit_fnb(x, y, weights, tau, result)`
Matches upstream `rq.wfit` scaling: rows of `x` and `y` are multiplied by
`weights` before fitting.

### `rq_fit_lasso(x, y, tau, lambda, result)`
`lambda` has length `p`. Set the intercept penalty to zero explicitly.

### `rq_fit_scad(x, y, tau, lambda, result)`
Iteratively reweighted SCAD QR corresponding to `rq.fit.scad`.

### `rq_fit_pfn(x, y, tau, result [, seed, ...])`
Dense Portnoy-Koenker preprocessing scheme for large `n`.

## Other algorithms

### `lprq(x, y, bandwidth, tau, ngrid, result)`
Local linear Gaussian-kernel QR.

### `nlrq_fit(y, theta0, tau, model, result)`
The callback is

```fortran
subroutine model(theta, fitted, jacobian)
```

and returns the fitted response and Jacobian at the current parameter vector.

### `kuantiles(x, probs, values [, qtype])`
Hyndman-Fan types 1 through 9; default is type 7.

### `qselect(x, prob)`
Order-statistic selection matching the package's `q489`/`qselect` convention
for ordinary interior probabilities.

### `rq_bootstrap_xy(x, y, tau, nrep, coefficients, info [, seed])`
Pairs bootstrap.

### `recursive_least_squares(...)`, `combinations(...)`, `random_exponential(...)`
Standalone translations of corresponding numerical helpers.
