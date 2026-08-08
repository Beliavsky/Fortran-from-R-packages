# API

## `subplex_control`

Fields:

- `reltol`: relative parameter tolerance; default `epsilon(1.0_dp)`.
- `maxeval`: maximum objective evaluations; default 10000.
- `parscale`: optional length-1 or length-`n` initial scale/step vector.

## `subplex_result`

Fields:

- `x`: estimated minimizer.
- `value`: objective value at `x`.
- `counts`: objective evaluations used by Subplex.
- `convergence`: `-2` invalid scale/input, `-1` evaluation limit, `0`
  tolerance satisfied, `1` machine-precision limit.
- `message`: diagnostic text corresponding to the convergence code.
- `hessian`: allocated when `compute_hessian=.true.`.

## `subplex_minimize`

```fortran
call subplex_minimize(fn, x0, result [, control] [, compute_hessian])
```

`fn` has the explicit interface

```fortran
function fn(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
end function fn
```

The objective is minimized. The R package does not expose maximization.

## `numerical_hessian`

```fortran
call numerical_hessian(fn, x, h, hess)
```

This reproduces the centered finite-difference Hessian used by the package's
C wrapper. `subplex_minimize` chooses `h = abs(parscale)*eps**(1/3)` when a
Hessian is requested.
