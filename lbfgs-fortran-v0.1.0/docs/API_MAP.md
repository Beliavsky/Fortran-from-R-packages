# R-to-Fortran API map

## Main routine

| R package | Fortran |
|---|---|
| `lbfgs(call_eval, call_grad, vars, ...)` | `call lbfgs_minimize(objective, gradient, x, result, parameters, ...)` |
| compiled combined evaluation internally | `call lbfgs_minimize(evaluate, x, result, parameters, ...)` |
| returned list `$value` | `result%value` |
| returned list `$par` | updated `x` argument |
| returned list `$convergence` | `result%status` |
| returned list `$message` | `result%message` |

## Parameter map

| R argument | `lbfgs_parameter_t` component |
|---|---|
| `m` | `m` |
| `epsilon` | `epsilon` |
| `past` | `past` |
| `delta` | `delta` |
| `max_iterations` | `max_iterations` |
| `linesearch_algorithm` | integer `linesearch` constant |
| `max_linesearch` | `max_linesearch` |
| `min_step` | `min_step` |
| `max_step` | `max_step` |
| `ftol` | `ftol` |
| `wolfe` | `wolfe` |
| `gtol` | `gtol` |
| machine epsilon `xtol` | `xtol` |
| `orthantwise_c` | `orthantwise_c` |
| `orthantwise_start` | `orthantwise_start`, converted to one-based indexing |
| `orthantwise_end` | `orthantwise_end`, inclusive in Fortran |

## Line-search constants

- `lbfgs_linesearch_morethuente`
- `lbfgs_linesearch_backtracking_armijo`
- `lbfgs_linesearch_backtracking_wolfe`
- `lbfgs_linesearch_backtracking_strong_wolfe`
- `lbfgs_linesearch_backtracking`, an alias for regular Wolfe and the required
  OWL-QN line search

## Callback interfaces

The combined evaluation callback computes objective and gradient together:

```fortran
subroutine evaluate(x, f, g, step, user_data)
```

The separate interface accepts:

```fortran
subroutine objective(x, f, user_data)
subroutine gradient(x, g, user_data)
```

An optional progress callback receives the current point, raw gradient,
objective, norms, accepted step, iteration number, and number of line-search
evaluations. Returning nonzero cancels the run.

## R-only infrastructure omitted

- R environments and `...` argument matching
- R function objects and Rcpp external pointers
- `.Call` registration
- R list construction
- console printing controlled by `invisible`

These are interface features rather than optimization algorithms.
