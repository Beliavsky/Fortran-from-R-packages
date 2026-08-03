# GARCHIto-fortran

A modern Fortran translation of the computational core of the R package
`GARCHIto` 0.1.0. The library estimates unified and realized GARCH-Ito
volatility models and returns fitted daily conditional variances together with
a one-step-ahead forecast.

The project is self-contained and uses no external numerical libraries.

## Implemented functionality

- `unified_est`: unified GARCH-Ito estimation from realized volatility and
  daily log returns
- `realized_est`: realized GARCH-Ito estimation with or without jump variation
- `realized_est_option`: realized GARCH-Ito estimation augmented by
  option-implied volatility, with homogeneous or power-heterogeneous errors
- Original parameter bounds and stationarity restriction
- Original quasi-likelihoods and unconditional initial variances
- OLS initialization for the option measurement equation
- Deterministic multistart bounded optimization
- Fitted conditional-variance series and one-step-ahead forecast

The upstream package has no package plotting function. Plotting shown only in
its vignette remains outside the numerical Fortran library. R package objects,
`.rda` loading, and vignette infrastructure are retained under `original/` but
are not translated.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example unified_fit
fpm run --example realized_fit
fpm run --example option_fit
```

The package uses standard Fortran 2018 and has no FPM dependencies.

## Minimal use

```fortran
use garchito, only : dp, garchito_result, realized_est

type(garchito_result) :: fit
real(dp) :: rv(100), jv(100)

! Fill rv and jv with nonnegative daily series.
call realized_est(rv, fit, jv=jv)
print *, fit%coefficients
print *, fit%pred
```

See [API.md](API.md) for the full interface and result fields.

## Optimizer and the Rsolnp dependency

The R package calls `Rsolnp::solnp`. This source distribution does **not** copy
or statically combine the earlier GPL-2-only Rsolnp Fortran translation,
because `GARCHIto` is GPL-3-only and those licenses are incompatible for a
single combined binary. Instead, this project contains an independently
implemented projected, bounded Nelder-Mead optimizer with deterministic
multistart and coordinate polishing. It enforces the same parameter bounds and
stationarity constraint used by the R package. See [PORTING.md](PORTING.md).

## License

The upstream package declares GPL-3. This translation is distributed under
**GPL-3.0-only**. Original metadata and computational R sources are retained
under `original/`.
