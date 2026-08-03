# API reference

## Types

### `nloptr_problem`

- `n`: number of parameters
- `n_ineq`, `n_eq`: numbers of inequality and equality constraints
- `lower`, `upper`: bound arrays
- `objective`: scalar callback
- `inequality`, `equality`: optional vector callbacks

### `nloptr_options`

Important fields include `algorithm`, `local_algorithm`, `ftol_rel`,
`ftol_abs`, `xtol_rel`, `xtol_abs`, `constraint_tol`, `maxeval`, `max_outer`,
`population`, `seed`, `initial_step`, `penalty_initial`, and `penalty_growth`.

### `nloptr_result`

Returns `solution`, `objective`, `max_constraint`, `status`, `iterations`,
`evaluations`, `converged`, and `message`.

## Main routine

```fortran
call nloptr(problem, x0, options, result)
```

The named wrappers have the same four main arguments, with a few optional flags:

```fortran
call lbfgs(problem, x0, options, result)
call varmetric(problem, x0, options, result, rank_one)
call tnewton(problem, x0, options, result, restart, precondition)
call auglag(problem, x0, options, result, derivative_free)
call stogo(problem, x0, options, result, randomized)
call mlsl(problem, x0, options, result, derivative_free)
```

## Numerical derivatives

```fortran
call nl_grad(x, objective, gradient, status)
call nl_jacobian(x, m, constraints, jacobian, status)
call check_derivatives(x, objective, analytic_gradient, result, tolerance)
```

## Options helpers

```fortran
options = nl_opts('NLOPT_LD_LBFGS')
options = nloptr_get_default_options()
call nloptr_print_options(options)
```
