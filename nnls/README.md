# nnls-fortran

Modern Fortran/FPM translation of the computational core of the R package
`nnls` 1.6.

The package solves

    min || A x - b ||_2

with either all coefficients constrained nonnegative (`nnls_fit`) or with a
mixture of nonnegative and nonpositive coefficients (`nnnpls_fit`).

## Features

- Lawson-Hanson active/passive-set NNLS algorithm.
- Mixed-sign NNNPLS constraints.
- Passive and bound index sets.
- Residuals, fitted values, residual norm and deviance.
- Lawson-Hanson-compatible status codes.
- Default `3*N` secondary-iteration limit.
- Free-form Fortran 2018 modules; no BLAS/LAPACK dependency.
- No implicit-interface callback machinery.

## Build

```text
fpm build
fpm test
```

The `scripts/test_gfortran.bat` and `scripts/test_gfortran.sh` scripts compile
with strict diagnostics independent of FPM.

## Basic use

```fortran
use nnls, only : dp, nnls_result, nnls_fit

type(nnls_result) :: fit
real(dp) :: a(m,n), b(m)

call nnls_fit(a, b, fit)
print *, fit%x
print *, fit%deviance
```

For mixed signs, `con(j) < 0` constrains `x(j) <= 0`; otherwise `x(j) >= 0`:

```fortran
call nnnpls_fit(a, b, con, fit)
```

See `API.md` and the programs under `example/`.
