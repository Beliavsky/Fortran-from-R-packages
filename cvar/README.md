# cvar-fortran

A modern Fortran 2018/FPM translation of the computational core of the R package **cvar 0.6**.

The library computes lower-tail value at risk (VaR) and expected shortfall (ES/CVaR) from:

- a quantile function;
- a cumulative distribution function;
- a probability density and matching quantile function;
- a sample vector; or
- each column of a sample matrix.

It also provides GARCH(1,1) simulation and forecasting with normal, standardized Student-t, and standardized generalized-error innovations.

## Main features

- Distribution-independent callback APIs
- Scalar and vector loss probabilities
- Empirical R type-7 quantiles
- Matrix-column VaR and ES
- Location-scale transformations
- Correct transformation from log returns to simple returns
- Adaptive Gauss-Kronrod tail integration
- Automatic CDF inversion by bracketing and bisection
- Normal, Student-t, standardized-t, and GED utilities
- Self-contained deterministic random-number generator
- GARCH(1,1) simulation with burn-in
- Analytical volatility forecasts
- Plug-in and Monte Carlo predictive intervals
- No external numerical-library dependency

## Build and test

```text
fpm build
fpm test
fpm run
fpm run --example garch_forecast
```

The project follows the standard FPM layout. It was also tested directly with GNU Fortran 14.2 because `fpm` was not installed in the translation environment.

## Small example

```fortran
program example
    use cvar, only : dp, normal_quantile, var_qf, es_qf
    implicit none

    print *, var_qf(normal_quantile, 0.05_dp)
    print *, es_qf(normal_quantile, 0.05_dp)
end program example
```

The expected values are approximately `1.64485362695` and `2.06271280751`.

## Callback convention

A custom distribution function has the scalar interface

```fortran
function distribution_function(x) result(value)
    use cvar, only : dp
    real(dp), intent(in) :: x
    real(dp) :: value
end function distribution_function
```

Module procedures are recommended. Host-associated internal procedures are legal Fortran but can require compiler-generated trampolines on some systems.

## Scope

Translated computational functionality:

- `VaR`, `VaR_qf`, and `VaR_cdf`
- `ES`
- numeric-vector and matrix methods
- `GarchModel`
- `sim_garch1c1`
- `predict.garch1c1`
- normal, standardized-t, and GED innovation support

Omitted R infrastructure includes S3 dispatch, argument-name matching, expression construction, R RNG-state preservation, package-loading logic, and R documentation machinery. Typed Fortran procedures replace those facilities.

See `PORTING.md` for the detailed mapping and intentional differences.

## License

`GPL-2.0-or-later`, preserving the original package's `GPL (>= 2)` declaration. See `LICENSE`, `NOTICE.md`, and `license/`.
