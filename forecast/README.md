# forecast-fortran

A modern free-format Fortran translation of the computational core of the R package
`forecast` 9.0.2, built for the Fortran Package Manager (FPM).

The upstream `forecast` package is GPL-3. This port preserves the upstream source tree
under `upstream/forecast-master/` and bundles the supplied Fortran translations of the
three computational dependencies used directly here:

- `fracdiff` for ARFIMA/fractional differencing;
- `urca` for unit-root based differencing selection;
- `nnet` for neural-network autoregression.

The combined package is distributed under GPL-3.0-only. See `NOTICE.md` and the
license files in each bundled dependency for provenance.

## Implemented numerical areas

- Box-Cox and inverse Box-Cox transformations, including bias adjustment and lambda selection.
- ACF, PACF, cross-correlation, tapered ACF/PACF with linear-process bootstrap intervals,
  Fourier terms, seasonal dummies, moving averages, differencing, and frequency/seasonal helpers.
- Forecast accuracy, Diebold-Marian testing, Ljung-Box diagnostics, and regression
  PRESS/CV/AIC/AICc/BIC/adjusted-R2 calculations.
- Mean, random-walk/naive, seasonal-naive, Croston/SBA/SBJ, and Theta forecasts.
- ETS state recursions translated from upstream `etscalc.c`, automatic ETS search,
  initial-state optimization, admissibility checks, fixed smoothing parameters, SES,
  Holt, Holt-Winters, analytic class-1/class-2/class-3 forecast variances, and
  simulation/bootstrap prediction intervals.
- ARIMA/SARIMA CSS initialization plus Gaussian ML fitting, diffuse integrated-state
  likelihoods, numeric ARIMAX/xreg, missing-observation ML filtering, fixed-structure
  refitting, simulation, impulse-response forecast variances, and expanded upstream-ordered
  `auto.arima`-style stepwise search including approximation truncation.
- `ndiffs` through the supplied `urca` implementation and OCSB/MSTL seasonal differencing.
- ARFIMA fitting/forecasting through the supplied `fracdiff` port.
- `nnetar`-style repeated neural-network autoregression through the supplied `nnet` port.
- BATS/TBATS state-space construction/filtering, seed-state estimation, parameter fitting,
  ARMA-error optimization, automatic model selection, forced control branches, structural
  refitting, non-integer TBATS seasonal periods, harmonic selection, and multi-step forecasts.
- Double-seasonal Holt-Winters (`dshw`).
- Cleveland-style LOESS/STL core and iterative MSTL decomposition, seasonal adjustment,
  trend/remainder extraction, and `stlf`-style reseasonalized forecasting.
- Matrix-based time-series linear regression (`tslm` numerical layer), linear `modelAR`,
  and stochastic cubic-spline forecasting with GCV smoothing selection.
- Missing-value interpolation, robust outlier detection/replacement, and time-series cleaning.
- Moving-block/decomposition bootstrap, bagged ETS, generic bagged-model callbacks, and
  blocked/unblocked `CVar`/rolling-origin CV callback APIs.
- Calendar helpers including month lengths, business-day counts with explicit holidays,
  and Gregorian Easter effects.

## Deliberately omitted R presentation/interface code

Formula parsing, S3 dispatch, `ts`/`msts` metadata, `ggplot2`/base plotting, graphical
geoms, model printing, data-frame wrangling, and bundled R datasets are not translated.
The Fortran APIs accept explicit vectors, matrices, periods, and model controls instead.

## Remaining parity differences

Version 0.3.0 closes the previously documented major numerical gaps. The remaining
differences are mostly implementation-equivalence or R-ecosystem concerns:

- integrated ARIMA uses a self-contained diffuse state-space implementation rather than
  R `stats::arima` internals, so likelihoods need not be bit-for-bit identical in difficult
  startup cases; ML supports missing observations, while automatic differencing/model
  selection still expects enough complete data for its unit-root tests;
- `auto_arima` follows the upstream stepwise candidate order and supports truncation and
  the main search controls, but does not reproduce R parallel execution or every edge-case
  warning/metadata branch;
- STL/MSTL and smoothing splines are self-contained Fortran implementations and are not
  bit-for-bit clones of R `stats::stl` / `smooth.spline`;
- named `timeDate` holiday-calendar datasets are not bundled; explicit holiday dates and
  Gregorian Easter effects are supported;
- R formula/S3/`ts` metadata, plotting, printing, and data-wrangling layers remain
  deliberately outside the numerical port.

See `PORTING_NOTES.md` and `API_MAPPING.md` for details.

## Build

```sh
fpm build
fpm test
fpm run --example demo_forecast
```

BLAS and LAPACK are required. The three supplied Fortran dependency ports are local FPM
path dependencies under `dependencies/`.
