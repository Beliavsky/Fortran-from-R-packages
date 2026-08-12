# API

All public symbols are re-exported by module `lpsolve`.

## Constants

- `dp`
- directions: `LP_MIN`, `LP_MAX`
- constraint senses: `LP_LE`, `LP_GE`, `LP_EQ`
- statuses: `LP_OPTIMAL`, `LP_SUBOPTIMAL`, `LP_INFEASIBLE`, `LP_UNBOUNDED`,
  `LP_NUMFAILURE`, `LP_TIMEOUT`

## Types

### `lp_control`

Important fields:

- `feasibility_tol`
- `optimality_tol`
- `integrality_tol`
- `max_simplex_iter`
- `max_nodes`
- `timeout_seconds`
- `scale_rows`
- `bland_rule`
- `num_binary_solutions`

### `lp_result`

- `status`
- `objective`
- `solution(:)`
- `duals(:)`
- `reduced_costs(:)`
- `solutions(:,:)` for requested multiple binary solutions
- `solution_objectives(:)`
- `solution_count`
- `simplex_iterations`
- `nodes`
- `sensitivity_ranges_available` (false in v0.1.0)

### `sparse_constraints`

Triplet sparse matrix with fields `nrow`, `ncol`, `row(:)`, `col(:)`, `val(:)`.

## General LP/MILP

```fortran
call solve_lp(direction, objective, A, sense, rhs, result [, control] &
              [, integer_variables] [, binary_variables])
```

`A` is `(n_constraints,n_variables)`.  Variables are nonnegative.  Integer and
binary variable arrays contain 1-based variable indices.

Sparse equivalent:

```fortran
call solve_lp_sparse(direction, objective, sparse, sense, rhs, result, ...)
```

## Assignment

```fortran
call lp_assign(cost, result [, direction] [, control] [, assignment])
```

If supplied, `assignment` is filled as a matrix with the same shape as `cost`.

## Transportation

```fortran
call lp_transport(cost, row_sense, row_rhs, col_sense, col_rhs, result, &
                  [direction] [, integer_variables] [, control] [, flow])
```

By default all transportation variables are integer, matching the R wrapper.

## 8 queens

```fortran
type(q8_triplets) :: q8
call make_q8(q8)
```

The result has 252 nonzeros in 42 constraints over 64 binary variables.

## Status text

```fortran
print *, trim(lp_status_message(result%status))
```
