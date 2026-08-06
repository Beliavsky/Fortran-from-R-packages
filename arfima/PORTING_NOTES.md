# Porting notes

## Independent replacement of `ltsa`

The R package imports GPL-licensed `ltsa` routines for Durbin-Levinson
likelihoods/residuals and Toeplitz forecasting. Because `arfima` itself is MIT,
this distribution does not copy or link the GPL implementation. The required
algorithms were independently implemented from their standard mathematical
recursions under the MIT license.

## Exact and equivalent computations

The following are direct translations of upstream formulas:

- AR/PACF transformations and admissibility ranges
- ARMA autocovariance linear equations
- FDWN, FGN, and PLA autocovariances
- integer/seasonal differencing and inverse differencing
- Box-Jenkins transfer-function recursion
- Durbin-Levinson exact concentrated likelihood
- covariance-based Gaussian simulation

The upstream seasonal/FARMA mixing uses an FFT circular convolution on long
truncations. This port uses the corresponding linear autocovariance convolution,
which avoids circular wraparound. Results agree away from the truncation
boundary but are not guaranteed bit-for-bit identical there.

## Fitting

R's `optim` calls are replaced by a self-contained Nelder-Mead optimizer.
AR/MA and transfer-denominator coefficients are optimized in unconstrained
PACF coordinates. Long-memory parameters use smooth logistic transformations:

- FD: `(-0.99, 0.49)`
- FGN Hurst: `(0.001, 0.999)`
- PLA alpha: `(0.001, 2.999)`

Covariance matrices use a finite-difference Hessian and a numerical Jacobian
back to natural parameter coordinates.

`fit_arfima_modes` provides multiple perturbed starts and likelihood sorting.
It does not reproduce R's exact Cartesian start-grid and optimizer-selection
heuristics.

The R `fixed` named-list convenience interface is not reproduced. Fixed models
can be evaluated directly with `lARFIMA`; fitting a partially fixed high-order
model requires constructing the reduced model explicitly.

## Information matrix

The upstream package contains specialized closed-form ARFIMA information
matrices. The Fortran implementation uses numerical derivatives of the log
spectrum obtained from the model autocovariance sequence. This also supplies a
consistent information calculation for FGN and PLA models, where the upstream
identifiability check is limited.

## Regression and transfer functions

Static regressors are differenced consistently with the response before
fitting time-series errors. Dynamic transfer functions support multiple input
series and the upstream concatenated `delta`/`omega` convention.

Simulation and exact forecasting from fitted dynamic transfer functions are
not implemented. The upstream `sim_from_fitted` function also rejects transfer
models. Static-regression simulation and forecasting accept replacement future
regressor matrices.

## R-only features omitted

- S3 classes and methods
- formulas, calls, data frames, `ts`, `xts`, and names
- plotting and graphical layouts
- parallel cluster management
- package datasets and package-news helpers
- formatted printing

The typed result structures retain all principal numerical outputs: parameters,
likelihood, variance, residuals, fitted values, AIC/BIC, covariance, standard
errors, and convergence diagnostics.
