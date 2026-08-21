# splines-fortran

A modern Fortran 2018 translation of the computational code in the R
package `splines` version 2.0-7 by Douglas M. Bates and William N. Venables.

The port provides dependency-free numerical routines for B-spline bases,
derivatives, regression-spline bases, natural and periodic interpolation,
piecewise-polynomial conversion, monotone inverse splines, and linear
interpolation. R formula handling, S3 dispatch, data frames, printing, and
plotting are intentionally omitted.

## Build with FPM

```sh
fpm test
fpm run --example basic_splines
```

## Build with GNU Fortran

On Unix-like systems:

```sh
./scripts/run_tests.sh
```

On Windows with `gfortran` in `PATH`:

```bat
scripts\run_tests.bat
```

## Main API

```fortran
use splines, only : dp, b_spline_t, poly_spline_t, spline_design, &
   bs_basis, natural_spline_basis, fit_interpolating_spline, &
   fit_periodic_spline, to_polynomial_spline, &
   inverse_monotone_spline, linear_interp
```

Spline objects expose scalar and vector evaluation through the type-bound
`evaluate` generic. An optional derivative order may be supplied.

## Example

```fortran
type(b_spline_t) :: curve
real(dp), allocatable :: fitted(:)
integer :: status

call fit_interpolating_spline(x, y, curve, status)
fitted = curve%evaluate(x_new)
```

See `example/basic_splines.f90` and `docs/API_MAP.md`.

## Licensing

The original package and this translation are distributed under GPL version
2 or, at your option, any later version. The original R and C computational
sources are retained under `original/` for attribution and traceability.
See `LICENSE` and `NOTICE`.

## Source provenance

Input archive: `splines-master.zip`

SHA-256: `7bc80ffe2f4f4f1496b505cdac31c5fc27b8e98beb4ea54c832b566f2cbe9f55`
