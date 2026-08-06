# ltsa-fortran

A modern Fortran 2018 translation of the computational code in the R package
`ltsa` 1.4.6.1 (Linear Time Series Analysis).

The library provides Durbin-Levinson recursions, Trench Toeplitz inversion and
likelihood calculations, exact linear forecasts, ARMA autocovariances,
Davies-Harte and recursion-based simulation, general linear-process filtering,
and innovation-variance estimation. R plotting, dynamic loading, and R object
infrastructure are not included.

## Build

With FPM:

```text
fpm test
fpm run --example ltsa_demo
```

With GNU Make:

```text
make check
make optimized
make example
```

The checked build uses Fortran 2018 conformance, warnings as errors, array and
bounds checks, floating-point traps, and NaN initialization.

## Minimal example

```fortran
program demo
    use iso_fortran_env, only : int64
    use ltsa
    implicit none

    type(ltsa_error) :: error
    type(forecast_result) :: forecast
    real(dp), allocatable :: r(:), z(:)
    real(dp) :: phi(1)

    phi = 0.8_dp
    call tacvf_arma(phi, [real(dp) ::], 104, 1.0_dp, r, error)
    call set_ltsa_seed(1234_int64)
    call dl_simulate(100, r, z, error)
    forecast = trench_forecast(z, r, 0.0_dp, 100, 5)

    print *, forecast%forecasts(1,:)
end program demo
```

## API style

Idiomatic names use underscores, such as `dl_acf_to_ar` and
`trench_loglikelihood`. The `ltsa_compat` module also supplies close R-name
wrappers with punctuation removed, such as `dlacftoar`, `trenchforecast`, and
`tacvfarma`.

All routines return either an `ltsa_error` or a typed result containing one.
No routine terminates the calling application for ordinary invalid input.

## License

GPL-2.0-or-later, matching the upstream `License: GPL (>= 2)` declaration.
See `LICENSE`, `licenses/GPL-2.0.txt`, and the preserved upstream source under
`upstream/`.
