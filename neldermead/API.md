# API

```fortran
use neldermead
```

## Objective callback

```fortran
function objective(x) result(f)
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function objective
```

## Unconstrained search

```fortran
call fminsearch(objective, x0, result)
```

To change defaults:

```fortran
type(nm_options) :: opt
opt = fminsearch_options(size(x0))
opt%max_iter = 1000
opt%tol_delta_f = 1.0e-8_dp
call fminsearch(objective, x0, result, opt)
```

## Bounds / Box complex

```fortran
call fminbnd(objective, x0, lower, upper, result)
```

Nonlinear inequality constraints use `c(x) >= 0`:

```fortran
subroutine constraints(x, c)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: c(:)
end subroutine constraints

call fminbnd(objective, x0, lower, upper, result, constraint=constraints, ncons=2)
```

## General search

```fortran
type(nm_options) :: opt
opt%method = 'fixed'       ! 'variable', 'fixed', or 'box'
opt%simplex0_method = 'spendley'
call neldermead_search(objective, x0, opt, result)
```

## Results

`nm_result` contains `x`, `f`, `iterations`, `func_count`, `restart_count`,
`status`, `converged`, final `simplex`, and optional best-point histories.
