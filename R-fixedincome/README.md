# R-fixedincome-fortran

A modern Fortran 2018 translation of the computational code in Wilson Freitas's
R package **fixedincome 0.0.5**. The project uses the Fortran Package Manager
(FPM), has no external numerical dependencies, and preserves the upstream MIT
license.

The library provides:

- simple, discrete, and continuous compounding and implied rates;
- terms and day-count conversion among days, months, and years;
- spot rates, spot-rate curves, forward rates, discount factors, and curve
  transformations;
- flat-forward, linear, log-linear, natural cubic, monotone Hermite, and Hyman
  style monotone interpolation;
- Nelson-Siegel and Nelson-Siegel-Svensson evaluation and bounded fitting;
- curve selection, exact/interpolated lookup, insertion, maturity dates, and
  interpolation error calculations.

Plotting, `ggplot2` integration, R S3/S4 display/indexing infrastructure,
bundled data, and external named holiday calendars are not translated.

## Build

```text
fpm build
fpm test
fpm run demo_fixedincome
```

Or use the included GNU Fortran scripts:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## Small example

```fortran
use fixedincome

type(spot_rate_curve_t) :: curve
real(dp), allocatable :: r(:)
integer :: status

curve = spotratecurve([0.049_dp, 0.051_dp, 0.053_dp], &
   term([30.0_dp, 180.0_dp, 365.0_dp], 'days'), &
   'discrete', 'actual/365', 'actual', status=status)
call set_interpolation(curve, interp_flatforward(), status)
r = interpolate(curve, [90.0_dp, 270.0_dp], status)
```

See `API.md`, `PORTING.md`, and the programs under `example/` for more detail.

## License

MIT. The original package metadata, R computational sources, tests, and
upstream license record are retained under `original/`.
