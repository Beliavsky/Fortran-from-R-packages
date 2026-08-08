# Porting notes

## Optimizer

The R package delegates global estimation to `rgenoud`. This port uses a
self-contained differential-evolution search followed by the bundled L-BFGS-B
3.0 implementation. This preserves bounded global/local likelihood estimation,
but random trajectories, evaluation counts, and convergence diagnostics do not
match `rgenoud` exactly.

## Kalman filter

The R fast path calls `FKF.SP`; the Fortran implementation is native and uses
Cholesky-based dense innovation solves. It accepts a different observation set
at every date and uses a Joseph covariance update for numerical stability.

A row with no observed contracts now propagates the predicted state and
covariance. The base-R loop in the supplied source left the state unchanged on
such rows.

## Simulation correction

The supplied `futures_price_simulate` constructs `H = ME^2` and then multiplies
a standard-normal draw by `H`, producing standard deviation `ME^2`. The Fortran
port multiplies by `ME`, consistent with the Kalman measurement covariance and
the package documentation.

## Forecast intervals

Futures forecast uncertainty is computed from the covariance accumulated over
the forecast horizon and the contract's maturity loadings. The supplied R
percentile branch instead uses covariance accumulated over the full contract
maturity, even when the forecast horizon is zero.

## American options

`LSMRealOptions::LSM_american_option` is replaced by a native
Longstaff-Schwartz implementation. Supported basis names are `power`,
`laguerre`, `hermite`, `legendre`, and `chebyshev`. The time grid values a
Bermudan approximation, converging toward an American value as `dt` decreases.

## Omitted R infrastructure

Printing, roxygen documentation machinery, data frames, row/column labels,
parallel R clusters, and R's S3/list conventions are not emulated. Plotting is
not part of the supplied package's exported computational API.
