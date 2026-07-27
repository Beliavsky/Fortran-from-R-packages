# garchsk-fortran

Modern Fortran/FPM implementation of the numerical algorithms in the R package
`GARCHSK` 0.1.0.

The library estimates GARCHSK and GJRSK models with time-varying conditional
mean, variance, skewness, and kurtosis. It is self-contained and has no external
numerical dependencies.

## Build

```text
fpm build
fpm test
fpm run garchsk_demo
fpm run --example basic_garchsk
fpm run --example estimate_and_forecast
```

A direct GNU Fortran validation script is provided for systems without FPM:

```text
sh scripts/validate.sh
```

On Windows with GNU Fortran:

```text
scripts\validate.bat
```

## Main API

```fortran
use garchsk

type(moment_path) :: path
type(forecast_result) :: fcst
type(estimate_result) :: fit

path = garchsk_construct(params, data)
value = garchsk_lik(params, data)
fit = garchsk_est(data)
fcst = garchsk_fcst(fit%params, data, 20)
```

The asymmetric model uses the corresponding `gjrsk_*` routines.

The public compatibility names mirror the R package:

- `skewness`, `kurtosis`
- `garchsk_construct`, `garchsk_lik`, `garchsk_ineqfun`, `garchsk_est`, `garchsk_fcst`
- `gjrsk_construct`, `gjrsk_lik`, `gjrsk_ineqfun`, `gjrsk_est`, `gjrsk_fcst`

## Model parameters

GARCHSK uses ten parameters:

```text
[a1, b0, b1, b2, c0, c1, c2, d0, d1, d2]
```

GJRSK uses thirteen parameters:

```text
[a1, b0, b1, b2, b3, c0, c1, c2, c3, d0, d1, d2, d3]
```

The `b2`, `c2`, and `d2` GJRSK terms multiply negative-shock indicators, while
the final coefficient in each block is the lagged conditional moment.

## Estimation

`garchsk_est` and `gjrsk_est` use a native constrained Nelder-Mead optimizer.
The returned `estimate_result` includes:

- parameters
- inverse-Hessian standard errors
- t statistics
- negative and ordinary log likelihoods
- corrected AIC and BIC
- upstream-compatible standard errors and information criteria
- convergence and evaluation diagnostics

The likelihood omits the constant `log(2*pi)` by default, matching the original
package. Set `include_constant=.true.` when estimating or evaluating a model to
include it.

## Forecasting

Forecast routines use the final observation and final conditional state for the
one-step forecast. Multi-step forecasts recursively propagate the model's
conditional moments.

## License and provenance

The project is licensed under `GPL-2.0-or-later`, matching the original package.
The unmodified R source and supplied archive are retained under `original/` and
`provenance/`. See `NOTICE`, `COVERAGE.md`, and `PORTING_NOTES.md`.
