# cccp-fortran

A modern Fortran/FPM translation of the computational algorithms in R package
`cccp` 0.3-3, *Cone Constrained Convex Problems*.

The library solves convex programs with:

- linear, quadratic, or callback-defined convex objectives;
- linear equality constraints;
- nonnegative-orthant constraints;
- second-order/Lorentz-cone constraints;
- positive-semidefinite matrix constraints; and
- callback-defined nonlinear convex inequalities.

It also implements the original package's geometric-programming, L1-regression,
and long-only risk-parity helpers.

## Build

With FPM:

```text
fpm test
fpm run
fpm run --example second_order_cone
```

The project links to BLAS and LAPACK. A reproducible GNU Fortran build is also
provided:

```text
./build_gfortran.sh debug
./build_gfortran.sh release
```

On Windows with GNU Fortran, link against an available BLAS/LAPACK build and use
an equivalent batch file or invoke FPM.

## Small example

```fortran
use cccp_api
real(dp) :: p(2,2), q(2), a(1,2), b(1)
type(cccp_solution) :: sol

p = 2.0_dp * reshape([2.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2,2])
q = [1.0_dp, 1.0_dp]
a = reshape([1.0_dp, 1.0_dp], [1,2])
b = 1.0_dp
sol = cccp(p, q, a, b)
print *, sol%x
```

## Scope

All computational R entry points are represented. Module `cccp_api` exposes the original generic name `cccp`; module `cccp` exposes the identical solver as `cccp_solve` to avoid a Fortran module-name collision. R reference classes, S4
method dispatch, formatted printing, and Rcpp module infrastructure are omitted.
The original source tree is retained under `original/` for provenance.

See `API.md`, `PORTING.md`, and `TESTING.md` for details.
