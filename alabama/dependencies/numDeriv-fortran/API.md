# API

Use the public umbrella module:

```fortran
use numderiv
```

The floating-point kind is `dp = kind(1.0d0)`.

## Callback interfaces

Real scalar callback:

```fortran
function f(x) result(value)
   real(dp), intent(in) :: x(:)
   real(dp) :: value
end function
```

Real vector callback:

```fortran
function f(x) result(value)
   real(dp), intent(in) :: x(:)
   real(dp), allocatable :: value(:)
end function
```

Complex-step variants use the same forms with `complex(dp)` arguments and
results. Complex-step functions must be analytic and must preserve complex
arithmetic internally.

## Options

```fortran
type(deriv_options) :: options
```

Fields:

- `eps`: additive step near zero; default `1e-4`
- `d`: relative initial step; default `1e-4`
- `zero_tol`: threshold selecting `eps`
- `r`: number of approximation levels; default `4`, with `r=1` supported
- `v`: step-reduction factor; default `2`
- `show_details`: print Richardson tables

`hessian_options()` returns the Hessian defaults, principally `d=0.1`.
`first_deriv_options()` returns the gradient/Jacobian defaults.

## grad

```fortran
call grad(func, x, gradient, method, side, options, status, message)
```

`method` is `"Richardson"` (default) or `"simple"`. `side` is an optional
integer vector with values `0`, `+1`, or `-1`; zero means a centered
Richardson derivative and means a forward step for the simple method.

For the R package's special elementwise-vector interpretation, use:

```fortran
call grad_elementwise(func, x, gradient, ...)
```

Complex-step forms are:

```fortran
call grad_complex(func, x, gradient, status, message)
call grad_elementwise_complex(func, x, gradient, status, message)
```

## jacobian

```fortran
call jacobian(func, x, jac, method, side, options, status, message)
call jacobian_complex(func, x, jac, status, message)
```

`jac` is allocated by the routine with shape `(size(func(x)), size(x))`.

## hessian

```fortran
call hessian(func, x, hess, options, status, message)
call hessian_complex(func, x, hess, options, status, message)
```

`hessian` uses the Bates-Watts/Richardson calculation. `hessian_complex`
calculates complex-step gradients and differentiates them with Richardson
extrapolation, matching the strategy of the R package. `hess` is allocated by
the routine.

## genD

```fortran
type(gend_result) :: result
call genD(func, x, result, options)
```

Fortran is case-insensitive, so `genD` and `gend` are the same name.
`result%dmat` has `m` rows and `p*(p+3)/2` columns. The first `p` columns are
the Jacobian. Remaining columns contain the lower triangle of each output
component's Hessian in the order `(1,1), (2,1), (2,2), (3,1), ...`.

Other fields are `f0`, `x`, `p`, `d`, `options`, `status`, and `message`.

## Status codes

- `nd_success`
- `nd_invalid_argument`
- `nd_nonfinite_value`
- `nd_shape_mismatch`
