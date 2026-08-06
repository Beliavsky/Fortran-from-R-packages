# tsgarch-fortran

A modern Fortran 2018 translation of the computational core of the R package
`tsgarch` 1.0.4.

The library implements univariate conditional-volatility filtering, estimation,
simulation, forecasting, inference, diagnostics, likelihood profiles, and
rolling value-at-risk backtests. R formulas, S3 classes, date-indexed containers,
parallel execution, tables, and plotting are intentionally not reproduced.

## Implemented volatility models

- GARCH
- GJR-GARCH
- APARCH
- EGARCH
- Family GARCH
- Component GARCH
- IGARCH
- EWMA

Arbitrary nonnegative ARCH and GARCH orders are supported where the model
permits them. Variance regressors, variance targeting, constant means, three
initialization methods, and model-specific stationarity constraints are also
available.

## Innovation distributions

The vendored `tsdistributions` translation supplies ten standardized laws:

- normal
- Student-t
- skew-normal
- skew-Student
- GED
- skew-GED
- NIG
- generalized hyperbolic
- Johnson SU
- generalized-hyperbolic skew-Student

## Main Fortran API

Most applications only need the umbrella module:

```fortran
use tsgarch
```

The central data types are:

- `type(garch_spec)`
- `type(garch_parameters)`
- `type(fit_options)`
- `type(garch_filter_result)`
- `type(garch_fit)`
- `type(garch_simulation)`
- `type(garch_forecast)`
- `type(backtest_result)`
- `type(profile_result)`

The principal procedures are:

- `validate_specification`
- `initialize_parameters`
- `filter_garch`
- `estimate_garch`
- `simulate_garch`
- `simulate_conditional`
- `forecast_garch`
- `persistence`
- `unconditional_variance`
- `news_impact`
- `probability_integral_transform`
- `covariance_opg`
- `covariance_sandwich`
- `confidence_intervals`
- `profile_likelihood`
- `backtest_var`

See `API_MAP.md` for the mapping from R functions and methods.

## Minimal example

```fortran
program example
  use ghyp_kinds, only : dp, i8
  use tsgarch
  implicit none

  real(dp) :: seed_data(30)
  type(garch_spec) :: spec
  type(garch_parameters) :: par
  type(garch_simulation) :: sim
  type(garch_fit) :: fit
  type(garch_forecast) :: forecast
  type(fit_options) :: options
  integer :: i

  seed_data = [(0.0_dp, i = 1, size(seed_data))]
  spec%model = 'gjrgarch'
  spec%distribution = 'std'
  spec%constant = .true.

  par = initialize_parameters(seed_data, spec)
  par%omega = 0.02_dp
  par%alpha = 0.05_dp
  par%gamma = 0.06_dp
  par%beta = 0.86_dp
  par%dist%shape = 8.0_dp

  sim = simulate_garch(spec, par, 300, burn=80, seed=20260804_i8)
  options%compute_inference = .false.
  fit = estimate_garch(sim%series(:,1), spec, options=options)
  forecast = forecast_garch(sim%series(:,1), fit, 5, paths=500, seed=17_i8)

  write(*,'(a,f12.4)') 'log likelihood: ', fit%log_likelihood
  write(*,'(a,5f10.5)') 'forecast sigma: ', forecast%sigma
end program example
```

## Building

With FPM:

```text
fpm test
fpm run --example demo_tsgarch
```

With GNU Make and GNU Fortran:

```text
make check
make optimized
make demo
```

The checked build uses Fortran 2018, all common warnings as errors, bounds and
runtime checking, and backtraces. The optimized build uses `-O3` and the same
strict language and warning settings.

## Optimizer and the supplied nloptr translation

The user-supplied `nloptr` Fortran translation is retained unchanged under
`provenance/dependencies/`, but is not compiled or linked. It is licensed
LGPL-3.0-or-later, while upstream `tsgarch` is GPL-2.0-only. Combining those
sources into one statically linked GPL-2-only library would create a license
conflict.

Estimation instead uses the bounded Nelder-Mead implementation already present
in the GPL-2-compatible `tsdistributions` translation. Equality and
stationarity restrictions are enforced through parameterization and objective
penalties. See `DEPENDENCY_NOTES.md` and `PORTING_NOTES.md`.

## License

The translated `tsgarch` and `tsdistributions` sources are distributed under
GPL-2.0-only. Vendored `ghyp` source files retain GPL-2.0-or-later notices. The
unlinked `nloptr` provenance archive retains its own LGPL-3.0-or-later license.

This is an independent language translation and is not affiliated with the
upstream package author.
