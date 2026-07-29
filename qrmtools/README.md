# qrmtools-fortran

Modern Fortran/FPM translation of the numerical core of the R package
`qrmtools` 0.0-19, *Tools for Quantitative Risk Management*.

The library is self-contained and uses no BLAS, LAPACK, optimization,
statistics, plotting, or R runtime dependency. Calculations use double
precision through

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Implemented numerical functionality

- GEV, GPD, Pareto, and GPD-tail density, distribution, quantile, and random
  functions.
- Truncated composite distributions assembled from normal, Student-t, GEV,
  GPD, or Pareto components.
- GEV quantile, probability-weighted-moment, and maximum-likelihood fitting.
- GPD method-of-moments, probability-weighted-moment, and maximum-likelihood
  fitting.
- Hill estimates, confidence intervals, sample and GPD mean-excess functions,
  and GPD tail-probability estimation.
- Nonparametric, Student-t, standardized Student-t, GPD, Pareto, and GPD-tail
  VaR and expected shortfall.
- Range VaR, geometric VaR, and geometric expectiles.
- Brownian motion, geometric Brownian motion, Brownian bridges, and inverse
  Brownian-increment recovery.
- Black-Scholes calls, puts, delta, gamma, vega, theta, rho, vanna, and vomma.
- Rearrangement and block-rearrangement algorithms for worst VaR, best VaR,
  and best ES.
- Fixed-grid and adaptive RA/ARA/ABRA-style bound construction through a
  margin-indexed quantile callback.
- Crude homogeneous bounds, Pareto/Wang bounds, and dual worst-VaR bounds.
- Elliptical and nonparametric Euler-style risk allocation.
- Mahalanobis and Mardia multivariate-normality diagnostics.
- Simple, logarithmic, and difference returns with exact inversion.
- Recursive hierarchical matrix construction.
- Reparameterized normal and standardized-Student GARCH(1,1) likelihood,
  native fitting, residuals, and innovation-based tail-index estimation.
- Native normal, Student-t, beta, gamma, and chi-square support routines;
  random generation; matrix algebra; numerical Hessians; and Nelder-Mead
  optimization.

## Basic use

```fortran
use qrmtools, only : dp, black_scholes, fit_gpd_mle, fit_result, &
  var_gpd, es_gpd

real(dp) :: excesses(8)
type(fit_result) :: fit

excesses = [0.15_dp, 0.24_dp, 0.38_dp, 0.51_dp, &
            0.79_dp, 1.10_dp, 1.55_dp, 2.30_dp]

fit = fit_gpd_mle(excesses, estimate_covariance=.false.)
if (.not. fit%ok) error stop trim(fit%message)

print *, black_scholes(0.0_dp, 100.0_dp, 0.05_dp, 0.20_dp, &
  100.0_dp, 1.0_dp)
print *, var_gpd(0.99_dp, fit%parameters(1), fit%parameters(2))
print *, es_gpd(0.99_dp, fit%parameters(1), fit%parameters(2))
```

## Building with FPM

```text
fpm build
fpm test
fpm run qrmtools_demo
fpm run --example evt_and_risk
fpm run --example simulation_and_bounds
```

The package uses standard automatic FPM discovery under `src`, `app`,
`example`, and `test`.

## Direct compiler validation

Linux and macOS with GNU Fortran:

```text
./scripts/validate.sh
```

Windows with GNU Fortran available on `PATH`:

```bat
scripts\validate.bat
```

## Typed result objects

Numerical procedures that can fail return typed results with an `ok` flag and
message rather than terminating the program. Important types include
`fit_result`, `hill_result`, `brownian_result`, `greeks_result`,
`rearrangement_result`, `ra_bounds_result`, `allocation_result`,
`test_result`, and `garch_result`.

## Scope exclusions

The numerical algorithms were translated. The following R-specific or
external-service infrastructure was not compiled:

- graphics and plot-building functions;
- financial-data downloads through `quantmod`, `xts`, and `zoo`;
- R condition-catching and list/S3 presentation helpers;
- the generic multi-series `fit_ARMA_GARCH` wrapper around external
  `rugarch` and `HoltWinters` objects;
- bundled data sets as compiled program inputs.

The self-contained reparameterized `fit_GARCH_11` algorithm is included.

## Threading note

The native Nelder-Mead wrappers for EVT, geometric risk measures, and GARCH
store their current objective data in module state. Separate fitting calls are
safe sequentially, but concurrent calls to the same fitter should be
externally serialized.

## Licensing and provenance

The upstream package is licensed under GNU GPL version 3 or later. This
translation preserves `GPL-3.0-or-later`, the original authorship, and the
complete GPLv3 text. The supplied source is retained unmodified under
`original/qrmtools-0.0-19`, and the uploaded archive and checksum manifests
are stored under `provenance`.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for detailed mapping
and validation information.
