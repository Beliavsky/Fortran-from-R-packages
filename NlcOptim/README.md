# NlcOptim modern Fortran

A modern Fortran 2018 translation of the computational code in the R package
`NlcOptim` 0.6. The library solves smooth nonlinear programs of the form

```text
minimize    f(x)
subject to  ceq(x) = 0
            c(x) <= 0
            Aeq x = Beq
            A x <= B
            lb <= x <= ub
```

The public `solnl` routine implements sequential quadratic programming (SQP)
with finite-difference or user-supplied derivatives, a damped BFGS Hessian,
an exact-penalty line search, and Goldfarb-Idnani quadratic subproblems from
the supplied `quadprog` Fortran translation.

## Features

- Nonlinear equality and inequality callbacks
- Linear equality and inequality constraints
- Lower and upper bounds
- Source-compatible forward finite differences by default
- Optional central finite differences
- Optional analytic objective gradients and nonlinear-constraint Jacobians
- Damped BFGS Hessian updates
- Elastic feasibility QP fallback when a linearized subproblem is inconsistent
- Structured Lagrange multipliers, status codes, diagnostics, and counters
- No BLAS, LAPACK, R, or C dependency

> Version 0.1.1 fixes an FPM duplicate-module collision between the application and example targets.

## Build

With FPM:

```sh
fpm test
fpm run
```

With GNU Make:

```sh
make check
make optimized
```

The checked target enables bounds, floating-point, and runtime checks. The
optimized target uses `-O3`.

## Minimal example

```fortran
program example
  use quadprog_kinds, only: dp
  use nlcoptim, only: nlc_result, solnl
  implicit none
  type(nlc_result) :: fit

  call solnl([8.0_dp, -7.0_dp], objective, fit)
  print *, fit%x, fit%objective
contains
  function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = (x(1) - 2.0_dp)**2 + 3.0_dp * (x(2) + 1.0_dp)**2
  end function objective
end program example
```

For code that must not request an executable stack, provide callbacks as
module procedures rather than internal procedures.

## License

The translated package is GPL-3.0-only, matching `NlcOptim`. The vendored
`quadprog` code remains GPL-2.0-or-later and is distributed here under GPL-3.
See `NOTICE.md` and the preserved upstream sources under `original/`.
