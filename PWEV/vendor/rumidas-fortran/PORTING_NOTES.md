# Porting notes

## Array orientation

Fortran MIDAS matrices have shape `(K+1, N)`, matching the mathematical and R
orientation.  Row 1 is the oldest included low-frequency observation and row
`K+1` is the contemporaneous observation.  The compatibility weight vector
ends in zero, so the contemporaneous value is not used.

## Dates and frequencies

The R routine `mv_into_mat` performs weekly, monthly, quarterly, or yearly date
matching with `xts` and `lubridate`.  The Fortran port deliberately separates
calendar handling from numerical code.  Supply, for every high-frequency
observation, the integer index of its current low-frequency period to
`lag_matrix_from_period_index` or `mv_into_mat`.

## Parameter ordering

Parameter vectors preserve upstream ordering.  Examples:

- non-skewed GM Normal: `alpha, beta, m, theta, w2`;
- skewed GM Student-t: `alpha, beta, gamma, m, theta, w2, shape`;
- skewed GMX Normal: `alpha, beta, gamma, z, m, theta, w2`;
- skewed DAGM Normal: `alpha, beta, gamma, m, theta_pos, w2_pos,
  theta_neg, w2_neg`.

The typed `*_parameter_count` functions should be used before allocating a
parameter vector for more complex two-MIDAS specifications.

## Likelihood compatibility

The Student-t expression in the upstream R source omits the constant
`-0.5*log(pi)`.  The Fortran translation intentionally retains that convention.
It does not affect maximum-likelihood parameter estimates, but absolute
log-likelihood, AIC, and BIC values differ from a fully normalized Student-t
log-density by a known constant per observation.

The Gaussian path follows the upstream use of the sample mean of returns.  The
Student-t path follows the upstream direct use of returns, which normally are
already demeaned financial returns.

## Optimization

`ugmfit` and `umemfit` use the attached `maxLik` Fortran package through an FPM
path dependency.  Bounds and the stationarity condition are explicit.  As in
the R package, multiple randomized candidates are evaluated and the best is
used as the optimizer start.

The likelihood callbacks keep a module-level copy of the active dataset because
Fortran procedure pointers do not carry arbitrary user data.  Simultaneous
fitting in multiple threads is therefore not supported.  Model evaluation after
fitting is reentrant.

Observation score matrices are computed by central finite differences and are
used for QML sandwich covariance estimates.  This is more expensive than the
basic Hessian covariance but preserves the intended upstream inference.

## Fitted quantities

For GARCH-MIDAS results:

- `conditional` is total conditional variance `g_t * tau_t`;
- `long_run` is long-run variance `tau_t`;
- `short_run` is the unit-mean short-run variance component `g_t`.

The low-level `*_cond_vol` compatibility routines return volatility, matching
their R names.  The R `ugmfit` routine squares those values before storing them;
the typed fit result does the same.

For MEM results, `conditional` is the positive conditional mean, `long_run` is
`tau_t`, and `short_run` is `mu_t`.

## Out-of-sample work

The R fitting wrappers automatically split `xts` objects.  In Fortran, callers
slice the in-sample arrays, fit once, and then evaluate the complete or held-out
arrays with the estimated coefficients.  This avoids hidden copies and date
semantics.

## Forecasts

`multi_step_ahead_pred` accepts the final short-run variance and final long-run
variance explicitly.  For X models it also accepts the final X value and an
AR(1) coefficient estimated by the caller.  This preserves the recurrence while
avoiding an embedded ARIMA implementation and R time-series objects.
