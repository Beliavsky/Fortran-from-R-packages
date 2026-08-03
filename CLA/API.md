# API

## Modules

### `cla`

Preferred descriptive API. It re-exports the public types, status constants,
and the following procedures:

- `critical_line(mu, covar, lower, upper [, tol_lambda, check_covariance])`
- `mean_sigma(weights, asset_mu, covar, sigma, mean_return)`
- `find_sigma(target_mu, result, covar [, equal_tolerance])`
- `find_mu(target_sigma, result, covar [, tolerance, equal_tolerance])`
- `mu_sigma_garch(prices [, arch_order, garch_order, distribution, max_iterations, tolerance])`

Distribution constants:

- `cla_distribution_normal`
- `cla_distribution_student`

### `cla_api`

Compatibility API corresponding to the exported R names:

- `CLA`: scalar-bound and vector-bound generic interfaces
- `MS`
- `findSig`
- `findMu`
- `muSigmaGarch`

## Result types

### `cla_result_t`

- `weights(n_assets,n_turning)`
- `free_mask(n_assets,n_turning)`
- `lambdas(n_turning)`
- `gammas(n_turning)`
- `sigma(n_turning)`
- `mu(n_turning)`
- `n_assets`, `n_turning`, `info`, `warnings`

### `cla_path_query_t`

- `value(:)`: sigma from `find_sigma`, or mu from `find_mu`
- `weights(:,:)`
- `info`

### `cla_garch_result_t`

- `mu(:)` and `covariance(:,:)`
- `forecast_sigma(:)`
- `conditional_sigma(:,:)`
- fitted `mean`, `omega`, `alpha`, `beta`, and Student-t `shape`
- optimizer convergence flags

## Status codes

- `cla_success = 0`
- `cla_invalid_input = 1`
- `cla_infeasible_bounds = 2`
- `cla_singular_system = 3`
- `cla_no_improvement = 4`
- `cla_out_of_range = 5`
- `cla_garch_failure = 6`
