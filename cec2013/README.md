# cec2013-fortran

Modern Fortran/FPM translation of the computational core of the R package
`cec2013` 0.1-5, which provides the 28 CEC-2013 real-parameter single-objective
benchmark functions.

## Features

- All 28 CEC-2013 benchmark functions.
- All upstream-supported dimensions: 2, 5, 10, 20, 30, 40, 50, 60, 70, 80,
  90, and 100.
- Official shift vectors and rotation matrices copied to `data/`.
- Reusable `cec2013_context` that loads a dimension's data only once.
- Scalar and batch evaluation.
- Published optimum-value helper.
- No R or external numerical-library dependency.
- Explicit modern Fortran interfaces and allocatable storage.

## Build

```sh
fpm build
fpm test
```

The project also contains strict compiler scripts:

```sh
./scripts/test_gfortran.sh
```

or on Windows:

```bat
scripts\test_gfortran.bat
```

## Basic use

```fortran
use cec2013
implicit none

type(cec2013_context) :: ctx
real(dp) :: x(10), f
integer :: status

x = 0.0_dp
call ctx%init(10, "data", status)
f = cec2013_evaluate(ctx, 1, x, status)
```

The data path is explicit because FPM does not standardize installation of
arbitrary runtime data files.  In this source distribution the benchmark data
are in `data/`.

For a matrix of points, `cec2013_evaluate_batch` uses Fortran layout
`x(n_dimension, n_points)`, so each **column** is a point.  This is the transpose
of the row-per-point convention used by the R wrapper.

## Fidelity

The translation intentionally preserves source-level behaviors that affect the
numerical benchmark, including some surprising expressions in the supplied C
implementation.  In particular:

- the Different Powers exponent uses C integer division;
- `asyfunc` leaves the destination element unchanged for nonpositive inputs;
- the supplied Griewank-Rosenbrock implementation computes a rotation and then
  overwrites it with the unrotated shifted vector;
- the first `10*n` values of `shift_data.txt` are read exactly as the C source
  does for dimension `n`.

These are not silently replaced by alternate textbook definitions.

See `TRANSLATION_COVERAGE.md` and `VALIDATION.md` for details.

## License

GPL-3.0-or-later, matching the upstream `GPL (>= 3)` declaration.  See
`LICENSE`, `NOTICE.md`, and `original/cec2013-master/`.
