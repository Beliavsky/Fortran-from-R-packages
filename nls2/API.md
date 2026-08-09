# API

## Main types

### `type(nls_control)`

Fields mirror the relevant `nls.control` concepts:

- `maxiter = 50`
- `tol = 1e-5`
- `min_factor = 1/1024`
- `scale_offset = 0`
- `central_diff = .false.`
- `warn_only = .false.`
- `diff_step = sqrt(epsilon)`

### `type(nls_result)`

Contains:

- `par`
- `linear_par` for partially-linear fits
- `fitted`
- weighted `residuals`
- `covariance`
- `rss`
- `sigma`
- `fin_tol`
- `iterations`
- `evaluations`
- `start_index`
- `status`
- `converged`

### `type(nls2_search_result)`

Contains:

- `best`
- `fits(:)` for all candidates
- `starts(:,:)`
- `start_rss(:)`
- `n_candidates`
- `status`

## Model callbacks

```fortran
subroutine nls_model(x, par, yhat, ierr)
subroutine nls_jacobian(x, par, jac, ierr)
subroutine plinear_basis(x, theta, basis, ierr)
```

For the partially-linear callback, predictions have the form

```text
yhat = basis(x, theta) * beta
```

where `theta` is nonlinear and `beta` is solved by weighted least squares at every evaluation.

## Main procedures

### `fit_nls`

```fortran
call fit_nls(model, x, y, start, result, control, weights, jacobian, lower, upper)
```

Gauss-Newton NLS with optional analytical Jacobian and bounds.

### `evaluate_model`

```fortran
call evaluate_model(model, x, y, par, result, weights)
```

Evaluate RSS at one parameter vector without optimization. This is the key operation used by `nls2` brute-force/random/LHS search modes.

### `fit_plinear`

```fortran
call fit_plinear(basis_fn, x, y, start_theta, n_linear, result, control, weights)
```

Variable-projection fit for models linear in `n_linear` coefficients conditional on nonlinear `theta`.

### `evaluate_plinear`

Evaluates a partially-linear model at a nonlinear parameter vector, solving only its linear coefficients.

### `nls2_fit`

```fortran
call nls2_fit(model, x, y, start_matrix, algorithm, result, control, ...)
```

Recognized algorithms:

- `default`
- `port` (bounded Gauss-Newton compatibility path)
- `brute-force`
- `grid-search`
- `random-search`
- `lhs`
- `CPoptim` (compatibility sampler; see coverage notes)

If `start_matrix` has two rows, those rows are interpreted as lower and upper limits for generated starts. Otherwise the rows are explicit starting values.

### `nls2_fit_plinear`

Recognized algorithms:

- `plinear`
- `plinear-brute-force`
- `plinear-brute`
- `plinear-random`
- `plinear-lhs`

### Search helpers

```fortran
call make_grid(lower, upper, maxiter, grid)
call random_uniform(lower, upper, points)
call latin_hypercube(lower, upper, points)
call seed_rng(seed)
```

`make_grid` intentionally follows the R source's `ceiling(maxiter^(1/k))` rule, so the number of grid points can exceed `maxiter`.

## Statistics helpers

```fortran
ll = nls_loglik(result, weights)
df = nls_df_residual(result, weights)
call pearson_residuals(result, r)
```

## Status codes

- `nls2_ok = 0`
- `nls2_maxiter = 1`
- `nls2_singular = 2`
- `nls2_bad_input = 3`
- `nls2_model_error = 4`
- `nls2_no_finite_start = 5`
