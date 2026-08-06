# FKF modern Fortran

A modern Fortran 2018 translation of the computational core of the R package
`FKF` 0.2.6. The library implements multivariate linear-Gaussian Kalman
filtering and Durbin-Koopman smoothing with constant or deterministic
 time-varying system matrices and arbitrary componentwise missing observations.

## Main API

```fortran
use fkf_module

type(fkf_model) :: model
type(fkf_result) :: filtered
type(fks_result) :: smoothed

call kalman_filter(model, y, filtered, corrected_missing_likelihood=.true.)
call kalman_smooth(model, y, filtered, smoothed)
```

Direct array-based compatibility entry points are also available:

```fortran
call fkf(a0, p0, dt, ct, tt, zt, hht, ggt, y, filtered)
call fks(model, y, filtered, smoothed)
```

Matrices may be constant (time extent 1) or time varying (time extent `n`).
Missing observations are represented by IEEE NaNs. For a partially observed
vector, only finite components enter the update. Output innovation, covariance,
inverse-covariance, and gain entries corresponding to missing components remain
NaN, matching the upstream package.

## Build

With FPM:

```text
fpm test
fpm run --example local_level_example
```

With GNU Make:

```text
make checked
make optimized
```

The checked build enables warnings as errors, bounds and runtime checks,
floating-point traps, and signaling-NaN initialization.

## License

The upstream package is GPL (>= 2). This translation is distributed under
`GPL-2.0-or-later`. See `licenses/`, `upstream/`, and `provenance/`.
