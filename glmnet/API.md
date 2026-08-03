# API

All public entities are re-exported by `use glmnet`.

## Types

### `glmnet_control_type`

Important fields:

- `alpha`: elastic-net mixing parameter in `[0,1]`
- `nlambda`: requested path length
- `lambda_min_ratio`: negative selects the glmnet-style automatic default
- `threshold`: convergence tolerance
- `max_iterations`: coordinate/proximal iterations
- `max_outer_iterations`: IRLS iterations
- `standardize`, `intercept`
- `grouped`: grouped multinomial penalty
- `probability_min`, `eta_max`

Use `default_glmnet_control()` and `update_glmnet_control()` when convenient.

### `glmnet_path_result`

A fitted path. The primary arrays are:

```text
lambda(nlambda)
intercept(nout, nlambda)
beta(nvars, nout, nlambda)
dev_ratio(nlambda)
objective(nlambda)
df(nlambda)
iterations(nlambda)
converged(nlambda)
```

`nout=1` for Gaussian, binomial, Poisson, and Cox models; it is the number of
responses or classes for multiresponse Gaussian and multinomial models.

### `glmnet_cv_result`

Contains the full-data fit, fold losses, cross-validation means and standard
errors, out-of-fold predictions, and `lambda_min`/`lambda_1se` selections.

### Other types

- `glmnet_assessment_result`
- `glmnet_roc_result`
- `glmnet_survival_data`
- `glmnet_sparse_csc`

## Fitting

### `fit_glmnet(x, y, family, result, ...)`

Generic interface for vector or matrix responses.

Vector response families:

- `'gaussian'`
- `'binomial'`
- `'poisson'`

Matrix response families:

- `'mgaussian'`
- `'multinomial'` (rows are class probabilities or counts)

Optional arguments are `control`, `weights`, `offset`, `lambda`,
`penalty_factor`, `lower`, `upper`, and `excluded`.

### Specialized fits

- `fit_multinomial_path(x, class_id, result, ...)`
- `fit_multinomial_matrix_path(x, y, result, ...)`
- `fit_cox_path(x, start, stop, event, result, ...)`
- `fit_custom_family_path(x, y, family_working, result, ...)`
- `fit_glmnet_sparse(x_csc, y, family, result, ...)`
- `big_glm(x, y, family, lambda, result, ...)`

`fit_cox_path` accepts optional strata and an `efron` flag. Otherwise it uses
Breslow ties.

A custom family callback has the interface:

```fortran
subroutine callback(y, eta, base_weight, working, irls_weight, deviance, status)
```

`gaussian_identity_working` is supplied as an example.

## Prediction and coefficients

- `predict_glmnet(fit, x, prediction, status, prediction_type, offset)`
- `predict_glmnet_at(fit, x, s, prediction, status, prediction_type, offset)`
- `coef_glmnet(fit, s, intercept, beta, status)`
- `nonzero_coef(fit, lambda_index, indices, status)`

Predictions have shape `nobs x nout x nlambda`. Supported prediction types are
`'link'`, `'response'`, `'class'`, and `'risk'` where applicable. Numeric
lambda interpolation is linear in log lambda, matching glmnet's usual path
interpolation convention.

## Cross-validation

- `cv_glmnet` for Gaussian, binomial, and Poisson models
- `cv_multinomial`
- `cv_cox`
- `build_predmat` copies the out-of-fold prediction array

Fold generation is deterministic when `seed` is supplied. `lambda_1se` is the
largest lambda within one standard error of the optimum.

## Relaxed fits

- `relax_glmnet` for Gaussian, binomial, and Poisson paths
- `relax_multinomial`
- `relax_mgaussian`
- `relax_cox`

Each active set is refitted with an effectively unpenalized ridge value of
`1e-10`, which stabilizes rank-deficient active designs.

## Assessment

- `assess_glmnet`
- `assess_multinomial`
- `assess_cox`
- `auc`
- `roc_glmnet`
- `confusion_glmnet`
- `glmnet_measures`
- `concordance_index`
- `cox_gradient`
- `coxnet_deviance`

## Data helpers

- `na_replace`
- `na_sparse_fix`
- `prepare_x`
- `make_x`
- `rmult`
- `stratify_surv`
- `dense_to_sparse`
- `sparse_to_dense`

Fortran cannot reproduce R data-frame factor expansion without a schema. The
helpers therefore operate on already numeric matrices.
