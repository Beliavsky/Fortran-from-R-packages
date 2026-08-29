# lbfgsb3-fortran

A modern Fortran 2018/FPM translation of the computational code in the R
package `lbfgsb3c`.

The package provides L-BFGS-B 3.0 for large bound-constrained optimization,
including the extra absolute/relative parameter stopping criteria exposed by
`lbfgsb3c`. The R, Rcpp, registration, and environment machinery has been
replaced by native Fortran procedure callbacks and typed result/control
objects.

## Main API

```fortran
use lbfgsb3_mod, only : dp, lbfgsb_control_t, lbfgsb_result_t, &
                       lbfgsb_minimize

real(dp) :: x(2), lower(2), upper(2)
type(lbfgsb_control_t) :: control
type(lbfgsb_result_t) :: result

x = [-1.2_dp, 1.0_dp]
lower = -2.0_dp
upper =  2.0_dp
control%pgtol = 1.0e-8_dp
control%reltol = 0.0_dp
call lbfgsb_minimize(x, objective_gradient, result, lower, upper, control)
```

The callback has the interface:

```fortran
subroutine objective_gradient(x, f, g, user_data)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: f
  real(dp), intent(out) :: g(:)
  class(*), intent(inout), optional :: user_data
end subroutine objective_gradient
```

`lbfgsb_minimize_fd` accepts an objective-only callback and computes a
bound-aware finite-difference gradient.

## Features

- L-BFGS-B 3.0 direct subspace-minimization algorithm
- lower-only, upper-only, two-sided, fixed, and unbounded parameters
- scalar or elementwise bound arrays
- `factr`, projected-gradient tolerance, memory size, and trace controls
- `lbfgsb3c` absolute and relative parameter-change stopping criteria
- separate function/gradient counts and raw reverse-communication diagnostics
- optional polymorphic user data
- optional monitor callback with cancellation
- non-finite objective/gradient detection
- objective-only bound-aware finite differences
- low-level `setulb` reverse-communication routine in `lbfgsb3_core_mod`

## Build

```text
fpm build
fpm test
fpm run --example bounded_rosenbrock
fpm run --example objective_only
```

GNU Fortran scripts are also provided in `scripts/`.

## Licenses

The translated wrapper and package are GPL-2.0-only. The L-BFGS-B numerical
kernel retains the supplied BSD 3-clause notice. See `NOTICE.md` and
`LICENSES/`.
