# Validation

Five deterministic test programs are included.

## `test_durbin`

- AR(1) conversion to final AR coefficients, PACF, and prediction variances
- Raw one-step residuals and innovation variances
- Concentrated likelihood against the closed AR(1) expression
- Recursion-based simulation with supplied innovations

## `test_toeplitz`

- Trench inverse against the closed AR(1) precision matrix
- Matrix-inverse product identity
- Bordering update against direct inversion
- Exact GLS mean for a constant series
- Equality of Trench and Durbin-Levinson concentrated likelihoods
- Exact-likelihood variance positivity

## `test_arma_forecast`

- AR(1) and MA(1) autocovariance fixtures
- Durbin-Levinson and exact prediction variances
- Updated and directly recomputed forecast tables
- Closed-form AR(1) forecasts

## `test_simulation`

- General linear-process convolution fixture
- Davies-Harte embedding condition
- Simulated AR(1) mean, variance, and lag-one covariance
- Source-compatible endpoint path
- Rejection of an invalid covariance sequence

## `test_innovation`

- AR/AIC innovation variance on a simulated AR(1)
- Kolmogoroff periodogram estimate
- Invalid short-series handling

Both `make check` and `make optimized` pass under GNU Fortran 14.2.0.
