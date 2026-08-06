# Source coverage

The translation targets the non-graphical computational behavior of `tsmarch`
1.0.3.

## Covered

### DCC (`R/dcc*.R`, `src/dcc.*`, `src/corfun.*`)

- model parameters and stationarity
- Gaussian and Student conditional likelihoods
- symmetric and asymmetric DCC recursions
- constant-correlation models
- filtering, estimation, simulation, and forecasting
- covariance/correlation extraction and numerical inference

### Copula-GARCH (`R/copula*.R`, `src/copula.*`)

- marginal probability transforms
- Gaussian and Student copula likelihoods
- dynamic and constant copula correlation
- estimation and simulation-based prediction

### ICA and GO-GARCH (`R/ica*.R`, `R/gogarch*.R`, `src/gogarch.*`,
`src/radical.*`)

- whitening
- FastICA and RADICAL-style ICA
- factor GARCH fitting
- covariance/correlation reconstruction
- coskewness and cokurtosis tensors
- simulation forecasts and portfolio moments

### Multivariate distributions and utilities (`R/mdistributions.R`,
`R/extra.R`, `R/utils.R`, `src/distributions.*`, `src/helpers.*`)

- multivariate normal and Student density/random generation
- EWMA and shrinkage covariance
- PSD repair and matrix conversions
- combinations and lagged matrices
- aggregation helpers

### Tests and risk (`R/tests.R`, computational portions of `R/methods.R`)

- Engle-Sheppard constant-correlation diagnostic
- VaR and expected shortfall from simulations or Gaussian moments
- portfolio path aggregation
- convolution-based density, CDF, and quantile evaluation

## Replaced by typed fields

The numerical content of `coef`, `fitted`, `residuals`, `logLik`, `vcov`,
`tscov`, `tscor`, `tscoskew`, and `tscokurt` is stored directly in `dcc_fit`,
`copula_fit`, `gogarch_fit`, and forecast result objects.

## Omitted

- `R/plots.R` and all graphics
- print, summary, and `flextable` formatting
- S3 method registration and R generic dispatch
- formula parsing, model frames, date/calendar handling, and `xts`/`zoo`
  metadata
- package datasets and documentation-generation infrastructure
- asynchronous/parallel R execution

These omissions do not remove standalone numerical model kernels.
