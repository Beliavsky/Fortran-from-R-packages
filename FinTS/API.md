# FinTS-fortran API

All public procedures and types are exported by:

```fortran
use fints
```

## Kinds and status codes

```fortran
integer, parameter :: dp = kind(1.0d0)
```

Status constants:

- `fints_ok`
- `fints_invalid_input`
- `fints_singular`
- `fints_nonstationary`
- `fints_iteration_limit`
- `fints_no_data`
- `fints_numerical_failure`
- `fints_io_error`

`fints_status_message(status)` returns a short allocatable character string.

## Summary statistics

```fortran
call FinTS_stats(x, result [, start])
```

Arguments:

- `x(:)`: observations; IEEE NaNs are omitted
- `start`: optional numeric start index or date code
- `result`: `summary_result`

The result contains `start`, `size`, `mean`, `standard_deviation`, `skewness`, `excess_kurtosis`, `minimum`, `maximum`, and `status`.

## ACF, covariance, and PACF

```fortran
call acf(x, result [, lag_max, acf_type, demean])
```

`acf_type` is one of `correlation`, `covariance`, or `partial`. Correlation and covariance results include lag zero. Partial-correlation results contain lags 1 through `lag_max`.

```fortran
call cross_acf(x, result [, lag_max, acf_type, demean])
```

For `x(observations, series)`, `result%value(lag+1,i,j)` contains the lagged covariance or correlation between series `i` and `j`.

## Ljung-Box and Box-Pierce tests

```fortran
call AutocorTest(x, result [, lag, test_type, degrees_freedom])
```

`test_type` is `Ljung-Box`, `Box-Pierce`, or `rank`. The result is a `test_result` with `statistic`, `degrees_freedom`, `p_value`, `lag`, `n_observations`, `method`, and `status`.

## ARCH LM test

```fortran
call ArchTest(x, result [, lags, demean])
```

The routine regresses squared observations on an intercept and their first `lags` lagged squares. The statistic is the regression R-squared times the number of fitted rows.

## Asymptotic PCA

```fortran
call apca(x, number_factors, result)
```

For `x(observations, series)`, the `apca_result` contains:

- `eigenvalues(:)`
- `factors(observations, number_factors)`
- `loadings(series, number_factors)`
- `r_squared(series)`
- `status`

The two-stage heteroscedastic scaling algorithm follows `R/apca.R`.

## ARMA theoretical ACF and roots

```fortran
call arma_true_acf(ar, ma, lag_max, result [, partial, complex_tolerance])
call plotArmaTrueacf(ar, ma, lag_max, result [, pacf, complex_eps])
```

The second name is a compatibility wrapper for the computational part of the original plotting routine. No plot is produced.

`arma_acf_result` contains:

- `roots(:)`: dynamic AR roots, inside the unit circle for a stationary model
- `value(:)`: ACF at lags 0 through `lag_max`, or PACF at lags 1 through `lag_max`
- `damping(:)` and `period(:)` for one representative of each complex-conjugate pair
- `stationary`, `partial`, and `status`

Low-level helpers:

```fortran
call findConjugates(x, representatives [, complex_tolerance])
call polynomial_roots(coefficients, roots, status [, tolerance, max_iterations])
```

Polynomial coefficients use ascending powers and a zero lower bound.

## ARIMA fitting

```fortran
call ARIMA(x, order, result [, seasonal_order, seasonal_period, xreg, &
   include_mean, transform_pars, method, n_cond, box_test_lag, &
   box_test_df, box_test_type, tolerance, max_iterations, initial])
```

Required arguments:

- `x(:)`: finite univariate observations
- `order(3)`: `[p,d,q]`
- `result`: `arima_result`

Optional arguments:

- `seasonal_order(3)`: `[P,D,Q]`
- `seasonal_period`: required to be at least 2 when a seasonal order is used
- `xreg(:,:)`: regressors with one row per observation
- `include_mean`: ignored when `d + D > 0`
- `transform_pars`: map unconstrained reflection coefficients to stationary/invertible AR/MA parameters
- `method`: `CSS`, `ML`, or `CSS-ML`; all use conditional Gaussian fitting
- `n_cond`: additional initial conditional observations to omit
- `box_test_lag`, `box_test_df`, and `box_test_type`
- `tolerance` and `max_iterations`
- `initial(:)`: initial optimizer-scale parameters in the order AR, MA, seasonal AR, seasonal MA, mean, regressors

`arima_result` includes fitted coefficients, conditional residuals, fitted values, variance, log likelihood, AIC, regression R-squared, convergence information, and a residual `box_test`.

## Interest and return conversion

```fortran
value = compoundInterest(interest [, periods, frequency, net_value])
log_return = simple2logReturns(simple_return)
```

Scalar and array forms are available. Infinite `frequency` requests continuous compounding. With `net_value=.false.` or omitted, `compoundInterest` returns gross value; with `.true.`, it returns the increase over one unit invested.

## Month conversion

```fortran
call as_yearmon2(x, result)
```

`x` may be a real array encoded as `yyyy.mm` or an integer array encoded as `yyyymm`. `yearmon_result` contains `year(:)`, `month(:)`, `converted`, `duplicate_count`, and `status`.
