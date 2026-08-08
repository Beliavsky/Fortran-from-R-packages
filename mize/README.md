# mize-fortran

A modern Fortran 2018 computational port of the R package `mize` 0.2.5.9000.
The package provides unconstrained gradient-based optimization through a native,
typed callback API and supports both one-shot and stateful iteration.

## Implemented algorithms

- Steepest descent (`SD`)
- Full-memory BFGS
- Symmetric rank-one (`SR1`)
- Limited-memory BFGS (`L-BFGS`)
- Nonlinear conjugate gradient with FR, CD, DY, HS, HS+, PR, PR+, LS, HZ,
  HZ+, and PR-FR updates
- Exact or finite-difference Newton directions
- Partial/frozen Hessian optimization (`pHess`)
- Truncated Newton / Newton-CG with exact Hessian-vector callbacks or finite
  differences
- Classical and Nesterov-style momentum
- Nesterov accelerated gradient schedules
- Delta-bar-delta per-coordinate learning rates
- Optional L-BFGS preconditioning for CG and truncated Newton
- Constant, Armijo backtracking, bold-driver, and safeguarded Wolfe searches
- Function- or gradient-based adaptive momentum restart
- Stateful `mize_init` / `mize_step` operation
- Convergence checks, progress storage, monitor cancellation, polymorphic user
  data, and finite-difference gradient diagnostics

## Build with FPM

```text
fpm build
fpm test
fpm run --example rosenbrock
fpm run --example stateful
```

## Minimal example

```fortran
program example
  use mize_mod, only : dp, mize_control_t, mize_result_t, mize_minimize
  implicit none

  type(mize_control_t) :: control
  type(mize_result_t) :: result
  real(dp) :: x(2)

  x = [-1.2_dp, 1.0_dp]
  control%method = 'L-BFGS'
  control%memory = 7

  call mize_minimize(x, rosenbrock_fg, result, control)
  print *, x, result%value

contains

  subroutine rosenbrock_fg(x, f, g, user_data)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: f, g(:)
    class(*), intent(inout), optional :: user_data

    f = 100.0_dp * (x(2) - x(1) ** 2) ** 2 + (1.0_dp - x(1)) ** 2
    g(1) = -400.0_dp * x(1) * (x(2) - x(1) ** 2) - 2.0_dp * (1.0_dp - x(1))
    g(2) = 200.0_dp * (x(2) - x(1) ** 2)
  end subroutine rosenbrock_fg
end program example
```

The objective callback returns the function and gradient together. Optional
Hessian, Hessian-vector, monitor, momentum-schedule, and user-data arguments are
available through `mize_minimize` and `mize_step`.

## Scope

This is a native computational port, not an emulation of R lists, closures,
method objects, or lifecycle hooks. The public numerical functionality is
represented by typed Fortran objects. See `docs/API_MAP.md` and
`docs/PORTING_NOTES.md` for exact mappings and implementation differences.

## License

BSD 2-Clause. The original source and attribution are retained in the release.
