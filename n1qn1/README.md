# n1qn1-fortran

Modern Fortran 2018 translation of the computational core of the R package
`n1qn1`, which ports Claude Lemarechal's Scilab `n1qn1` unconstrained
full-memory BFGS optimizer.

The translation preserves the original packed LDL-transpose Hessian update and
Wolfe-type line search. It does not substitute a generic BFGS implementation.

## Features

- Full-memory BFGS optimization for smooth unconstrained objectives
- Original `n1qn1` line-search logic
- Analytic objective and gradient callbacks
- Optional polymorphic user data
- Optional progress/cancellation callback
- User-supplied initial Hessian
- Exact restart from an internal packed LDL-transpose factor
- Returned dense Hessian, packed full Hessian, and internal factor
- Original workspace-size helper and R-compatible padded `c_hess`
- Explicit status codes and non-finite-value validation
- No external numerical dependencies

Because the method stores a dense symmetric Hessian approximation, memory and
work per iteration grow quadratically with the number of parameters. It is
intended for small and medium-sized smooth problems rather than large-scale
optimization.

## FPM build

```text
fpm build
fpm test
fpm run --example banana
```

## Basic use

```fortran
use n1qn1_module, only : dp, n1qn1_result_t, n1qn1_minimize

type(n1qn1_result_t) :: result
real(dp) :: x0(3)

x0 = [1.02_dp, 1.02_dp, 1.02_dp]
call n1qn1_minimize(objective, gradient, x0, result)
```

Callbacks have explicit interfaces:

```fortran
function objective(x, user_data) result(value)
    use n1qn1_module, only : dp
    real(dp), intent(in) :: x(:)
    class(*), intent(inout), optional :: user_data
    real(dp) :: value
end function objective

subroutine gradient(x, g, user_data)
    use n1qn1_module, only : dp
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    class(*), intent(inout), optional :: user_data
end subroutine gradient
```

Module procedures are recommended for callbacks.

## Reusing curvature

To provide a conventional dense Hessian estimate:

```fortran
call n1qn1_minimize(objective, gradient, x0, result, &
                    initial_hessian=previous%hessian)
```

To restart from the exact internal LDL-transpose representation:

```fortran
call n1qn1_minimize(objective, gradient, x0, result, &
                    initial_factor=previous%factor)
```

These two forms are deliberately distinct. See `docs/PORTING_NOTES.md`.

## Result fields

- `x`: final parameter vector
- `value`: final objective value
- `gradient`: final gradient
- `hessian`: final dense Hessian approximation
- `factor`: packed internal LDL-transpose representation
- `c_hess`: packed dense Hessian followed by zero workspace padding, matching
  the length used by the R package
- `iterations`: iterations attempted
- `function_evaluations`, `gradient_evaluations`: callback counts
- `gradient_norm_squared`: squared Euclidean norm of the final gradient
- `status`, `message`: termination information

## Validation

The tests reproduce the supplied R package's three-variable Rosenbrock case,
including:

- 40 objective/gradient evaluations with automatic Hessian initialization
- The reference final dense Hessian to approximately `2e-9`
- 29 evaluations with the documented user-supplied Hessian
- 33 evaluations when reusing the first fit's dense Hessian

Additional tests cover a positive-definite quadratic, a one-variable problem,
exact factor restart, evaluation limits, user data, and callback cancellation.

## Scope

Rcpp evaluation wrappers, R environments, external pointers, dynamic
registration, R printing controls, and assignment into R environments are not
translated. Their numerical role is replaced by native Fortran callbacks and
typed result objects.

## License

CeCILL version 2.1. See `LICENCE`, `COPYING`, and `NOTICE.md`.
