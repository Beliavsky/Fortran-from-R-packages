# TIMSAC Modern Fortran

This project is a standalone modern Fortran port of the numerical core in the
R package `timsac` 1.3.8-6.

The original package contains about 29,600 lines of fixed-form Fortran and
about 3,800 lines of R interface, allocation, printing, and plotting code. This
port converts the complete numerical core to free-form Fortran and provides a
small modern API for the most commonly used operations.

## Modernization performed

- fixed-form `.f` source converted to free-form `.f90`
- source identifiers and executable code normalized to lower case
- `double precision` declarations replaced by `real(dp)`
- `dp` defined as `kind(1.0d0)`
- `implicit none` added to every procedure
- fixed-form continuations converted to free-form continuations
- spaced legacy tokens such as `E N D` and `. LE.` normalized
- all converted files compile with `gfortran -std=f2018`
- an fpm project, Makefile, tests, and an example program are included

The numerical algorithms still intentionally retain some historical control
flow and storage constructs, including labeled `do` loops, `goto`, and
`common` blocks. They are valid to gfortran in Fortran 2018 mode but are not a
complete idiomatic redesign of every internal routine.

## High-level API

Use module `timsac` for these allocation-safe interfaces:

- `autocorrelation`
- `multivariate_correlation`
- `power_spectrum`
- `matrix_filter`
- `white_noise`

`matrix_filter` is a direct Fortran translation of the R-only `mfilter`
function. The other functions wrap the converted TIMSAC numerical routines.

Example:

```fortran
program demo
  use timsac, only: dp, autocorrelation, autocorrelation_result
  implicit none

  real(dp) :: x(5)
  type(autocorrelation_result) :: result

  x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
  result = autocorrelation(x, max_lag=2)

  print *, result%mean
  print *, result%covariance
  print *, result%correlation
end program demo
```

## Complete converted numerical core

All original compiled entry routines and their support procedures are in
`src/legacy`. The file `API_MAP.md` maps every exported R function to its
converted Fortran entry point. The modern module `timsac_raw` currently gives
explicit interfaces to the routines used by the high-level API. Other entry
points can be wrapped similarly using their declarations in `src/legacy`.

## Build with Make

```text
make
make test
make example
```

The library is created as `build/libtimsac.a`.

For a checked debug build:

```text
make debug
```

## Build with fpm

```text
fpm build
fpm test
fpm run
```

## Directory layout

```text
src/timsac_kinds.f90   floating-point kind
src/timsac_raw.f90     explicit interfaces used by the modern API
src/timsac.f90         modern allocation-safe API
src/legacy/*.f90       complete converted TIMSAC numerical core
app/                   runnable example
test/                  regression tests
tools/                 fixed-form conversion utility
reference/             original package metadata and copyrights
```

## License and attribution

The source is derived from the TIMSAC R package, which declares GPL version 2
or later. This modern Fortran port is distributed under GPL-2.0-or-later; the
complete license text is in `LICENSE`. Original copyright and package metadata
are retained in `reference/`. See `NOTICE` and `ORIGIN.md` for authorship,
modification, and provenance details.
