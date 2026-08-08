# API

## Optimizers

```fortran
call trust_optim(x0, fn, gr, "SR1", result [, control])
call trust_optim(x0, fn, gr, "BFGS", result [, control])
call trust_optim(x0, fn, gr, hs, result [, control])

call trust_optim_sr1(x0, fn, gr, result [, control])
call trust_optim_bfgs(x0, fn, gr, result [, control])
call trust_optim_sparse(x0, fn, gr, hs, result [, control])
```

`fn` is a scalar function with interface

```fortran
function fn(x) result(f)
  real(dp), intent(in) :: x(:)
  real(dp) :: f
end function
```

`gr` is

```fortran
subroutine gr(x, g)
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: g(:)
end subroutine
```

The sparse Hessian callback is

```fortran
subroutine hs(x, h)
  real(dp), intent(in) :: x(:)
  type(sparse_symmetric_matrix), intent(inout) :: h
end subroutine
```

`h` stores the lower triangle only.  `h%set_from_dense(A)` is provided for
small problems and testing; large applications should fill `row`, `col`, and
`val` directly in lower-triangle coordinate format.

## Controls

`type(trustoptim_control)` mirrors the R controls:

- `start_trust_radius`
- `stop_trust_radius`
- `cg_tol`
- `prec`
- `report_freq`
- `report_level`
- `report_precision`
- `report_header_freq`
- `maxit`
- `contract_factor`
- `expand_factor`
- `contract_threshold`
- `expand_threshold_ap`
- `expand_threshold_radius`
- `function_scale_factor`
- `precond_refresh_freq`
- `preconditioner`
- `trust_iter`

`function_scale_factor < 0` performs maximization exactly as in the R package.
SR1 forces an identity preconditioner when preconditioner 1 is requested,
matching the R wrapper.

## Result

`type(trustoptim_result)` contains

- `solution`
- `fval` (unscaled user objective)
- `gradient` (unscaled user gradient)
- `hessian` for sparse mode
- `iterations`
- `status`
- `trust_radius`
- `nnz`
- `method`
- `hessian_update_method`
- last CG iteration count/reason

Use `result%status_message()` and `result%method_name()` for readable labels.

## Sparse matrix

`type(sparse_symmetric_matrix)` provides

- `clear()`
- `matvec(x,y)`
- `to_dense(A)`
- `set_from_dense(A [, drop_tol])`

Only the lower triangle is stored.  Off-diagonal entries are reflected during
matrix-vector multiplication and dense conversion.

## Binary-choice example

Module `trustoptim_binary` provides

```fortran
type(binary_data)
type(binary_priors)
binary_value(...)
binary_gradient(...)
binary_hessian(...)
```

Both upstream parameter orderings (`order.row = FALSE/TRUE`) are supported.
