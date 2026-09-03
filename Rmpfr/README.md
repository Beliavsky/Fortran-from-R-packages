# Rmpfr — computational-core translation to modern Fortran

This directory is a modern free-form Fortran translation of the portable
computational core of the R package **Rmpfr 1.1-2**.  Rmpfr is fundamentally an
interface to the GNU MPFR arbitrary-precision floating-point library, which in
turn uses GMP.  This translation therefore binds directly to MPFR with
`ISO_C_BINDING` instead of replacing arbitrary precision by `real(dp)`.

The public `rmpfr` module provides an owning `mpfr_real` scalar type, deep-copy
assignment, per-value binary precision, directed rounding, MPFR exponent-range
controls, arithmetic/comparison operators, elementary and special functions,
constants, exact/integer combinatorics, probability kernels, summaries,
multiple-precision matrix products, Romberg integration, one-dimensional root
finding and optimization, qnorm inversion, and Hooke-Jeeves optimization.

## External numerical requirements

No BLAS, LAPACK, or ARPACK library is used.  The only native libraries are the
same essential system dependencies as upstream Rmpfr:

- GNU MPFR >= 3.2.0
- GNU GMP >= 4.2.3

They are declared as FPM native link dependencies in `fpm.toml`.

All ordinary Fortran real variables use `dp = real64`.  The binding module uses
`real(c_double)` only where the MPFR C ABI itself requires a C `double`.

### Windows with gfortran

A convenient supported-style setup is the MSYS2 UCRT64 toolchain.  In an
UCRT64 shell, MPFR and GMP can be installed with:

```text
pacman -S mingw-w64-ucrt-x86_64-mpfr mingw-w64-ucrt-x86_64-gmp
```

Use the matching UCRT64 gfortran/FPM toolchain, and ensure the UCRT64 `bin`
directory remains on `PATH` when running executables so the MPFR/GMP DLLs can
be located.  No separately installed BLAS or LAPACK is required.

## Build and test

From the extracted top-level `Rmpfr` directory:

```text
fpm build
fpm test
fpm run --example basic_usage
```

The current validation environment did not contain `fpm`; see
`BUILD_VALIDATION.md` for the strict direct-gfortran validation that was run
instead and for the literal FPM command attempts.

## Minimal use

```fortran
use rmpfr

type(mpfr_real) :: x, y

x = mpfr_from_string('0.1', 256)
y = mpfr_exp(x)
print *, trim(mpfr_to_string(y, 70))
```

`mpfr_real` values own MPFR storage and assignments are deep copies.  Matrix
and cumulative routines use output subroutines rather than allocatable array
function results; this avoids compiler-dependent ownership problems for arrays
of finalizable derived types.

## Scope

This is a computational translation, not an R compatibility layer.  R S4/S3
classes and method dispatch, R vector recycling, R printing/formatting UI,
serialization of R objects, and conversion to/from the R `gmp` package's
`bigz`/`bigq` objects are intentionally not reproduced.  See
`API_COVERAGE.md` for detailed coverage and intentional differences.

## Licensing and provenance

The translation is distributed under GPL-2.0-or-later, matching upstream
Rmpfr.  MPFR and GMP are external libraries and are not vendored.  Upstream
metadata is retained under `upstream/`; see `NOTICE.md` for attribution and
algorithm provenance.
