# lsei-fortran

Modern Fortran/FPM translation of the computational core of the R package
`lsei` 1.3-1.

The package solves least-squares and quadratic-programming problems with
nonnegativity, equality, inequality and box constraints.  It is standalone and
requires no R, BLAS or LAPACK installation.

## Main API

```fortran
use lsei

type(ls_result) :: fit
call nnls_solve(a, b, fit)
call pnnls_solve(a, b, kfree, fit)
call ldp_solve(e, f, fit)
call lsi_solve(a, b, e, f, fit)
call lsei_solve(a, b, c=c, d=d, e=e, f=f, res=fit)
call qp_solve(q, p, fit, c=c, d=d, e=e, f=f)
call pnnqp_solve(q, p, kfree, fit)
```

`hfti_solve`, `indx`, and `mat_maxs` are also provided.

## Build

```text
fpm build
fpm test
```

For a compiler-independent strict check on Windows, run
`scripts\\test_gfortran.bat`.

## Licensing

The R package declares `GPL (>= 2)`.  The bundled Lawson-Hanson source code is
identified by the upstream package as public domain.  The complete upstream
source tree is retained under `original/lsei-master/`.
