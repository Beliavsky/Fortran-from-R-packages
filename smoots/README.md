# smoots-fortran

A modern Fortran 2018/FPM translation of the computational algorithms in the R package **smoots 1.1.4**.

The original package estimates smooth deterministic trends and their derivatives in equidistant time series by local polynomial regression with iterative plug-in bandwidth selection. It also provides kernel smoothing, confidence bounds, ARMA forecasting, bootstrap intervals, and rolling forecast evaluation.

## Scope

The port includes the numerical content of every exported computational function:

- `gsmooth`: fixed-bandwidth local polynomial trend and derivative estimation
- `knsmooth`: Nadaraya-Watson kernel smoothing
- `tsmooth`: configurable iterative plug-in bandwidth selection
- `msmooth`: all ten predefined algorithms `A`, `B`, `N`, `NA`, `O`, `OA`, `OM`, `NM`, `OAM`, and `NAM`
- `dsmooth`: data-driven first- and second-derivative estimation
- `confBounds`: asymptotically unbiased estimates and pointwise normal confidence bounds
- `rescale`: derivative rescaling from `[0,1]` to actual equidistant time units
- `critMatrix` and `optOrd`: ARMA information-criterion matrices and constrained order selection
- `normCast`: Gaussian ARMA forecast intervals
- `bootCast`: forward-bootstrap ARMA forecast intervals
- `trendCast`: linear or constant trend extrapolation
- `modelCast`: combined nonparametric-trend and ARMA-residual forecasts
- `rollCast`: holdout forecasting, interval breaches, MASE, and RMSSE
- `maInfty`: infinite-MA coefficient recursion
- Bühlmann's iterative lag-window long-run variance estimator
- AR-, MA-, and ARMA-based long-run variance estimators
- low-level point forecast, simulated true forecast, covariance, and smoothing kernels formerly implemented in Rcpp

Plotting, colors, progress bars, parallel-worker management, S3 printing, and R `ts`/data-frame attributes are presentation or runtime facilities and are not included.

## Build with FPM

```console
fpm build
fpm test
fpm run
fpm run --example fixed_smoothers
fpm run --example derivative_example
```

No external Fortran library is required.

## Minimal example

```fortran
program example
   use smoots
   implicit none

   integer, parameter :: n = 240
   real(dp) :: y(n), x
   integer :: i
   type(smooth_result) :: fit

   do i = 1, n
      x = real(i, dp) / real(n, dp)
      y(i) = 2.0_dp + x + 0.3_dp * sin(12.0_dp * x)
   end do

   call msmooth(y, fit, p=1, mu=1, algorithm='A')
   print '(a,f10.6)', 'bandwidth = ', fit%b0
   print '(a,f10.6)', 'last trend = ', fit%estimate(n)
end program example
```

## Main API conventions

R lists are represented by typed derived types:

- `smooth_result`
- `arma_model`
- `confidence_result`
- `forecast_result`
- `rolling_result`

Status codes are returned in each result object and are also available as the `sm_ok`, `sm_invalid_input`, `sm_singular`, `sm_iteration_limit`, and `sm_fit_failed` constants.

## Numerical compatibility

The local polynomial and kernel smoothing recursions, boundary weighting, hidden lookup tables, plug-in exponents, variance-bandwidth multipliers, derivative kernels, lag-window recursion, MA-infinity recursion, and forecasting formulas follow the original source.

The R package delegates ARMA estimation to `stats::arima`, which uses R's exact/conditional likelihood and Kalman machinery. This self-contained port instead uses conditional Gaussian least squares with an analytical residual Jacobian and Levenberg-Marquardt iterations. Consequently, algorithms using `Mcf='AR'`, `'MA'`, or `'ARMA'`, automatic ARMA order selection, and bootstrap re-estimation can differ slightly from R. See `PORTING.md`.

## License

The original package declares `GPL-3`. This translation is distributed under **GPL-3.0-only**. Original source and metadata are retained under `original/` for attribution and provenance.
