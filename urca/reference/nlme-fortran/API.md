# API

Import the public API with:

```fortran
use nlme
```

## Kinds and status

- `dp = kind(1.0d0)`
- `NLME_SUCCESS`
- `NLME_INVALID_ARGUMENT`
- `NLME_SINGULAR`
- `NLME_NOT_POSITIVE_DEFINITE`
- `NLME_MAX_ITER`
- `NLME_NONFINITE`
- `NLME_CALLBACK_ERROR`
- `NLME_DIMENSION_ERROR`
- `nlme_status_message(status)`

## Model specifications

### `correlation_spec`

Fields: `kind`, `p`, `q`, allocatable `par`, `nugget`, and `fixed`.

Kinds:

- `COR_NONE`
- `COR_AR1`
- `COR_CAR1`
- `COR_ARMA`
- `COR_COMPOUND_SYMM`
- `COR_EXPONENTIAL`
- `COR_GAUSSIAN`
- `COR_LINEAR`
- `COR_RATIO`
- `COR_SPHERICAL`
- `COR_UNSTRUCTURED`

`par` uses natural parameters in the public interface. Examples are an AR(1)
correlation in `(-1,1)`, a CAR(1) base correlation in `(0,1)`, AR then MA
coefficients for `COR_ARMA`, and range followed by optional non-nugget ratio for
spatial correlations.

### `variance_spec`

Kinds: `VAR_CONSTANT`, `VAR_FIXED`, `VAR_IDENT`, `VAR_POWER`,
`VAR_EXPONENTIAL`, `VAR_CONST_POWER`, and `VAR_CONST_PROP`.

### `pd_spec`

Kinds: `PD_IDENT`, `PD_DIAG`, `PD_LOG_CHOL`, and `PD_COMPOUND_SYMM`.
The parameters are unconstrained log-scale/Cholesky parameters.

### `nlme_control`

Controls iteration limits, convergence tolerance, finite-difference step,
REML/ML selection, covariance optimization, and tracing.

## Core fitting procedures

### `fit_gls`

```fortran
call fit_gls(y, x, result, correlation, variance, time, group, &
             var_covariate, var_group, coordinates, control)
```

Fits a Gaussian GLS model. Correlation and variance parameters are estimated
unless their `fixed` component is true or `control%optimize_covariance` is false.

### `fit_lme`

```fortran
call fit_lme(y, x, z, group, result, random, correlation, variance, &
             time, var_covariate, var_group, coordinates, control, fixed_sigma)
```

Fits a Gaussian linear mixed model from explicit design matrices. The random
covariance, residual scale, and residual covariance parameters are estimated by
ML or REML. `fixed_sigma` reproduces the fixed-residual-scale use case.

### `fit_gnls`

```fortran
call fit_gnls(model, y, xdata, theta0, result, correlation, variance, ...)
```

The callback interface is:

```fortran
subroutine model(theta, x, yhat, status)
  real(dp), intent(in) :: theta(:), x(:,:)
  real(dp), intent(out) :: yhat(:)
  integer, intent(out) :: status
end subroutine
```

### `fit_nlme`

```fortran
call fit_nlme(model, y, xdata, group, theta0, random_index, result, ...)
```

`random_index` gives the nonlinear parameter positions that receive group-level
random effects.

## Result types

- `gls_result`
- `lme_result`
- `nonlinear_result`

They contain parameter estimates, covariance matrices, fitted values,
residuals, random effects where applicable, scale, log likelihood, diagnostics,
iteration count, status, and convergence flag.

## Covariance helpers

- `correlation_matrix`
- `arma_autocorrelation`
- `spatial_distance_matrix`
- `variance_sd`
- `pd_matrix`
- `build_residual_covariance`
- `gls_log_likelihood`

## Numerical and diagnostic helpers

- `finite_difference_gradient`
- `finite_difference_hessian`
- `numerical_model_jacobian`
- `acf_values`
- `empirical_variogram`
- `pooled_sd`
- `simulate_lme`
- `is_balanced`
- `group_summary`
- `fit_lm_list`
- `fit_nls_list`
