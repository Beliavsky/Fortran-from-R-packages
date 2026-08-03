# API reference

All public procedures are available through:

```fortran
use fincovregularization
```

Real calculations use `real(dp)`. Procedures that can fail accept an optional
integer `status`; `fincov_ok` means success. `fincov_status_message(status)`
returns a readable description.

## Norms

### `f_norm2(matrix)`

Returns the squared Frobenius norm, `sum(matrix**2)`.

### `o_norm2(matrix [, status])`

Returns the squared spectral/operator norm. It is computed as the largest
eigenvalue of `transpose(matrix) * matrix`.

## Covariance regularization

### `banding(sigma [, k, status])`

Keeps entries with `abs(i-j) <= k`; the default is `k=0`.

### `tapering(sigma, l [, h, status])`

Applies the tapering weights used by the R package. The default ratio is
`h=0.5`.

### `hard_thresholding(sigma [, threshold, status])`

Sets off-diagonal entries with absolute value below the threshold to zero. The
default threshold is `0.5`.

### `soft_thresholding(sigma [, threshold, status])`

Shrinks off-diagonal entries toward zero while preserving the diagonal. The
default threshold is `0.5`.

### `ind_cov(sigma [, status])`

Matches the original `Ind.Cov`: computes the sample covariance of the rows of
`sigma` as observations and retains only its diagonal.

### `threshold_min(sigma [, method, tolerance, status])`

Finds the R package's minimum hard- or soft-threshold constant by bisection on
the smallest eigenvalue. `method` defaults to `"hard"`.

## Cross-validation

The three selectors return a `type(cv_result)` with fields:

- `regularization`, `method`
- `parameter_opt`, `parameter_index`
- `parameter_grid(:)`, `cv_error(:)`
- `n_cv`, `norm`, `seed`, `h`, `status`

The split is repeated deterministically with a portable Park-Miller random
number generator. The training size follows the R package:
`ceiling(n * (1 - 1/log(n)))`.

### `banding_cv(data [, n_cv, norm, seed])`

Searches integer bandwidths from zero through `p-1`.

### `tapering_cv(data [, h, n_cv, norm, seed])`

Searches taper lengths from zero through `p-1`.

### `threshold_cv(data [, method, thresh_len, n_cv, norm, seed])`

Searches an evenly spaced hard- or soft-threshold grid. `norm` accepts `"F"`
or `"O"` case-insensitively.

## Factor covariance models

Input return matrices use shape `(n_observations, n_assets)`.

### `macro_factor_cov(assets, factor_or_factors [, status])`

Generic interface accepting either a vector or matrix of macro factors. Fits
multivariate OLS with an intercept and returns factor covariance plus diagonal
residual variance.

### `fundamental_factor_cov(assets, exposure [, method, status])`

`exposure` has shape `(n_assets, n_factors)`. `method` is `"OLS"` or `"WLS"`
and defaults to WLS.

### `stat_factor_cov(assets [, k, status, selected_k])`

Uses an eigen decomposition equivalent to the original centered-data SVD. With
`k=0` or omitted, the original singular-value-greater-than-one rule selects the
number of factors. `selected_k` optionally reports it.

## Portfolio routines

### `gmvp(covariance [, allow_short, status, tolerance, max_iterations])`

Returns global minimum-variance weights summing to one. Short sales are allowed
by default. The long-only case uses a convex active-set method and enforces
nonnegative weights.

Unlike the R display-oriented function, the Fortran routine does not round
weights to four decimal places.

### `risk_parity(covariance [, status, tolerance, max_iterations, objective_value])`

Minimizes the original package's equal-risk-contribution objective with
Nelder-Mead. As in the R routine, the last weight is `1-sum(first weights)` and
no explicit nonnegativity constraint is imposed.

### `risk_parity_objective(weights, covariance)`

Evaluates the translated objective directly.
