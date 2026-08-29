# Porting notes

## Design

The port uses explicit vector/matrix APIs rather than attempting to emulate R's formula,
S3, `ts`, `msts`, list, or data-frame object systems. The public umbrella module is
`forecast`; specialized modules can also be used directly.

All newly written Fortran is free form and uses `implicit none`. `dp` is defined from
`real64`. BLAS/LAPACK are used for least-squares and dense linear algebra.

## Dependency handling

The user supplied Fortran translations of `fracdiff`, `urca`, and `nnet`. They are bundled
unchanged under `dependencies/` and referenced by FPM path dependencies. The forecast
translation calls them only where upstream `forecast` actually needs the corresponding
numerics:

- `fracdiff`: ARFIMA fit, fractional differencing and fractional weights.
- `urca`: KPSS test used by `ndiffs`.
- `nnet`: neural networks used by `nnetar`.

## Native BATS/TBATS translation

The observation, innovation and transition matrices follow the indexing of the upstream
Rcpp/Armadillo sources. During validation an off-by-one error in the first translation of
the AR/MA observation/transition blocks was identified by comparing the zero-based C++
indices with Fortran one-based indices; it has been corrected.

High-level BATS/TBATS seed states are estimated from the upstream relation
`D = F - g w^T`. A zero-seed filter produces the transformed residual sequence and the
successive `w D^(t-1)` rows form the seed-state regression. BATS seasonal seed states use
the sum-to-zero constraint; TBATS trigonometric states do not need that constraint.

## ARIMA

The upstream R package delegates core ARIMA fitting to R's `stats::arima`/`stats::arima0`.
Version 0.3.0 uses CSS as an initializer/optional approximation and a self-contained
Gaussian likelihood for ML fitting. Integrated models use explicit inverse-differencing
states with diffuse initialization; missing observations skip the measurement update while
the state propagates. Numeric xreg coefficients are fitted jointly with ARMA parameters,
future regressors are required for ARIMAX forecasts, and `arima_refit` reuses a fitted
structure/coefficients on new data. Regular and seasonal AR/MA polynomials are combined
multiplicatively with an explicit coefficient origin, the process-mean parameterization
follows the R convention, and forecast uncertainty uses ARMA plus differencing impulse
responses. `auto_arima` follows the upstream candidate/restart order and supports
approximation truncation with a full-sample ML refit. Bit-for-bit equivalence to the R
`stats::arima` initialization is not claimed.

## ETS

`ets_state_forecast`, `ets_update`, and `ets_calc` are translations of the upstream native
state recursions in `etscalc.c`. Version 0.3.0 optimizes initial states, enforces the
upstream-style forbidden-model and admissibility constraints, supports fixed smoothing
parameters correctly, expands automatic model search, implements class-1/class-2/class-3
analytic forecast variances, and adds normal/bootstrap simulation intervals using R type-8
quantiles.

## Spline

The stochastic cubic-spline covariance from upstream `spline.R` is retained. GCV fitted
values now treat the constant-plus-linear spline null space exactly through a universal
smoother projection rather than a finite large-prior-variance approximation. Forecast
covariances retain the finite construction used by upstream `forecast::spline_model`. The
implementation remains self-contained and is not bit-for-bit `stats::smooth.spline`.

## Cleaning/decomposition

Version 0.2.0 includes a self-contained Cleveland-style LOESS/STL core with robust bisquare
iterations and iterative MSTL extraction. It is used by seasonal-strength differencing and
STLF-style reseasonalization. Because it does not call R's `stats::stl`, bit-for-bit parity
with R's implementation is not claimed.

## BATS/TBATS refitting

Fitted BATS/TBATS objects retain their seed state and training series. `bats_refit` and
`tbats_refit` hold the fitted state-space structure and coefficients fixed, rerun the
filter on new observations, and update fitted values, residuals, terminal state, and
innovation variance. Automatic fitting also exposes force/exclude/auto choices for Box-Cox,
trend, and damping controls.

## Missing observations

Stationary and integrated ARIMA ML likelihoods accept IEEE NaN observations. Missing
measurements skip the Kalman update while state prediction continues. CSS fitting remains
restricted to complete observations, and automatic unit-root/seasonal-differencing tests
are not advertised as a general missing-data imputation layer.
