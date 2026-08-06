# quadprog-fortran

Modern Fortran/FPM translation of the computational code in R package
`quadprog` 1.5-8.

The library solves strictly convex quadratic programs of the form

```text
minimize    0.5*x^T*D*x - d^T*x
subject to  A(:,1:meq)^T*x  = b(1:meq)
            A(:,meq+1:q)^T*x >= b(meq+1:q)
```

The implementation preserves the Goldfarb-Idnani dual active-set algorithm
used by the upstream package. Both the dense and compact indexed-column
solvers were modernized from fixed-form Fortran into free-form Fortran 2018
modules.

## Features

- Dense `solve_qp` interface corresponding to R `solve.QP`
- Compact `solve_qp_compact` interface corresponding to `solve.QP.compact`
- Equality and inequality constraints
- Optional pre-factorized Hessian input as the upper-triangular `R^{-1}` for
  `D = R^T R`
- Returned unconstrained solution, objective, Lagrange multipliers, active
  constraints, and iteration counters
- Explicit status codes instead of R exceptions or process termination
- No BLAS, LAPACK, or other external dependency

## Build with FPM

```sh
fpm build
fpm test
fpm run --example basic_qp
fpm run --example portfolio_qp
fpm run demo_quadprog
```

The package name is `quadprog-fortran`; the Fortran module is `quadprog`.

## Minimal use

```fortran
use quadprog_kinds, only: dp
use quadprog, only: qp_result, solve_qp

type(qp_result) :: fit
fit = solve_qp(dmat, dvec, amat, bvec, meq=1)
if (.not. fit%succeeded()) error stop fit%message
print *, fit%solution
```

See `API.md`, `PORTING.md`, and the `example/` directory for details.

## License

GPL-2.0-or-later, matching the upstream package's `GPL (>= 2)` license.
The original package sources and license files are retained under `original/`.
