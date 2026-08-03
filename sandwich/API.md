# API

Import the complete public interface with `use sandwich`.

All real calculations use `dp = kind(1.0d0)`. Allocatable outputs are allocated
by the called procedure. Optional integer status arguments return one of:

- `SANDWICH_SUCCESS`
- `SANDWICH_INVALID_ARGUMENT`
- `SANDWICH_DIMENSION_MISMATCH`
- `SANDWICH_SINGULAR_MATRIX`
- `SANDWICH_INSUFFICIENT_DATA`
- `SANDWICH_NUMERICAL_FAILURE`
- `SANDWICH_UNSUPPORTED`

## Generic sandwich calculations

### `meat(scores, meat_matrix, status, adjust)`

Computes `transpose(scores) * scores / n`. With `adjust=.true.`, multiplies by
`n/(n-k)`.

### `sandwich_covariance(bread, meat_matrix, nobs, covariance, status)`

Computes `bread * meat * bread / nobs` and symmetrizes the result.

### `vcov_opg(scores, covariance, status, adjust)`

Returns the inverse of `transpose(scores) * scores`, optionally multiplied by
`n/(n-k)`.

### `bread_from_information(information, nobs, bread, status)`

Returns `nobs * inverse(information)`.

## Linear regression

### `type(ols_model)`

Contains `nobs`, `ncoef`, `df_residual`, `sigma2`, `coefficients`, `fitted`,
`residuals`, `weights`, `bread`, classical `covariance`, observation `scores`,
and `hat` values.

### `fit_ols(x, y, model, status, weights, offset)`

Fits ordinary or weighted least squares. `x` is `n x k`. Weights must be
nonnegative. The offset is subtracted before fitting and added to fitted values.

### `ols_scores`, `ols_bread`, `ols_hatvalues`

Standalone OLS adapter procedures.

## Heteroscedasticity-consistent covariance

### `hc_weights(residuals, hat, df_residual, type, omega, status)`

Accepted types are `const`, `HC`, `HC0`, `HC1`, `HC2`, `HC3`, `HC4`, `HC4m`,
and `HC5`.

### `meat_hc(x, residuals, type, meat_matrix, status, hat, omega)`

Computes the HC meat matrix. A caller may provide leverage values or completely
custom nonnegative variance weights `omega`.

### `vcov_hc(x, residuals, bread, type, covariance, status, hat, omega)`

Computes the complete robust covariance.

## HAC covariance and kernels

### `kernel_weight` and `kernel_weights`

Supported kernels are truncated/uniform, Bartlett/triangular, Parzen,
Tukey-Hanning, and quadratic spectral/QS.

### `prewhite_var(scores, order, residuals, recolor, status)`

Fits a zero-intercept VAR of the requested order to the score matrix and returns
prewhitened residuals and the recoloring matrix.

### `meat_hac` and `vcov_hac`

Use an explicit lag-weight vector. Element 1 is the lag-zero weight. Optional
`adjust` defaults to true and `prewhite_order` defaults to zero.

### `bandwidth_andrews` and `andrews_weights`

Andrews automatic bandwidth selection using an AR(1) or conditional-sum-of-
squares ARMA(1,1) approximation. `prewhite_order` defaults to one.

### `bandwidth_newey_west` and `newey_west_weights`

Newey-West automatic bandwidth and Bartlett lag weights. Automatic bandwidth
supports Bartlett, Parzen, and quadratic-spectral kernels.

### `lumley_weights`

Constructs truncate or smooth WEAVE weights from the residual isotonic ACF.

### `long_run_variance`

Computes the robust covariance of a sample mean using Andrews or Newey-West
weights, corresponding to the upstream `lrvar()` calculation.

### `type(hac_diagnostics)`

Contains bias correction, approximate degrees of freedom, bandwidth, effective
sample size, number of retained weights, and prewhitening order.

## Clustered covariance

### `meat_cluster(scores, cluster, meat_matrix, status, type, cadjust, multi0, x, residuals, hat)`

`cluster` is an `n x p` integer matrix. Every nonempty subset of the `p`
clustering dimensions is combined by inclusion-exclusion. Types `HC0`, `HC1`,
`HC2`, and `HC3` are supported. `cadjust` defaults to true. `multi0` substitutes
the ordinary HC0 component for the highest-order intersection. Cluster HC2 and
HC3 require the OLS design matrix and residuals.

### `vcov_cluster`

Wraps `meat_cluster` and the generic sandwich calculation. `fix=.true.` projects
an indefinite multiway result to the positive-semidefinite cone.

## Panel covariance

### `meat_panel_longitudinal` and `vcov_panel_longitudinal`

Panel longitudinal/Driscoll-Kraay covariance. Inputs include an observation
score matrix, integer cluster IDs, and integer time IDs. Available lag constants
are `PANEL_LAG_MAX`, `PANEL_LAG_NW1987`, and `PANEL_LAG_NW1994`. With
`aggregate=.true.` scores are summed across clusters by time; false retains
within-cluster serial products.

### `meat_panel_corrected` and `vcov_panel_corrected`

Panel-corrected covariance from an OLS design matrix, residuals, cluster IDs,
and time IDs. `pairwise=.false.` uses complete time periods; true estimates each
cross-sectional covariance from its available pairwise periods.

## Bootstrap and jackknife

### `bootstrap_covariance(replicates, covariance, status)`

Sample covariance of a matrix of replicate estimates.

### `jackknife_covariance(estimates, covariance, status, center_estimate)`

Returns `(R-1)/R` times the cross-product of deviations. The default center is
the replicate mean; a full-sample estimate may be supplied.

### `vcov_bootstrap_ols`

Cluster-resampling covariance for OLS. Supported types are `xy`, `jackknife`,
`fractional`, `residual`, `rademacher`, `mammen`, `norm`, and `webb`. Multiple
cluster dimensions are combined by inclusion-exclusion. `seed` gives repeatable
Fortran RNG initialization and `fix=.true.` projects the result to the PSD cone.

## Auxiliary numerical routines

The umbrella module also exports PAVA/isotonic ACF procedures and the matrix
routines used internally: linear solves, inversion, covariance, symmetric
Jacobi eigendecomposition, symmetric matrix powers, and PSD projection.
