# API

All public functionality is re-exported by module `rsolnp`.

## Callback interfaces

```fortran
subroutine objective_callback(x, value, data)
subroutine gradient_callback(x, gradient, data)
subroutine vector_callback(x, value, data)
subroutine jacobian_callback(x, jacobian, data)
```

`x` is a rank-one `real(dp)` parameter vector. The optional polymorphic `data`
argument receives `problem%data` when that component is allocated. Constraint
Jacobians use shape `(number_of_constraints, number_of_parameters)`.

## `type(solnp_problem)`

Important components:

- `name`: descriptive problem name.
- `n`: number of parameters; may be left zero and inferred from `start`.
- `n_eq`, `n_ineq`: numbers of equality and inequality functions.
- `fn`, `gr`: objective and optional analytic gradient callbacks.
- `eq_fn`, `eq_jac`: equality function and optional Jacobian callbacks.
- `ineq_fn`, `ineq_jac`: inequality function and optional Jacobian callbacks.
- `data`: optional polymorphic callback data.
- `start`, `lower`, `upper`: initial parameters and parameter bounds.
- `eq_b`: equality targets, so the solver enforces `eq_fn(x) = eq_b`.
- `ineq_lower`, `ineq_upper`: two-sided bounds on `ineq_fn(x)`.
- `best_fn`, `best_par`: optional benchmark reference solution.

Call `prepare_problem` to validate and fill omitted default bounds/targets.

## `type(solnp_control)`

- `rho`: initial augmented-Lagrangian penalty, default `1`.
- `max_iter`: maximum major iterations, default `400`.
- `min_iter`: maximum projected-BFGS inner iterations, default `400`.
- `delta`: finite-difference step scale, default `1e-6`.
- `tol`: convergence tolerance, default `1e-8`.
- `trace`: zero for quiet operation; positive values print iteration progress.
- `penalty_growth`: penalty multiplication factor, default `5`.
- `max_rho`: maximum penalty, default `1e12`.
- `line_search_max`: maximum Armijo backtracking iterations.
- `armijo`: Armijo sufficient-decrease coefficient.
- `min_step`: minimum line-search step.
- `restoration_iter`: feasibility-restoration iteration limit.

## `type(solnp_result)`

- `pars`, `objective`: final parameters and original objective value.
- `objective_history`, `constraint_history`, `step_history`: major-iteration histories.
- `lagrange`: estimated equality and active-inequality multipliers.
- `hessian`: final BFGS approximation to the augmented-objective Hessian.
- `ineq_slack`: final bounded slack variables.
- `out_iterations`, `inner_iterations`, `n_eval`, `elapsed`.
- `convergence`, `message`.
- `kkt`: `kkt_diagnostics` result.

Status constants are `solnp_success`, `solnp_max_iterations`,
`solnp_invalid_problem`, `solnp_numerical_failure`, and `solnp_infeasible`.

## Single-start solvers

```fortran
call solnp(problem, result [, control])
call csolnp(problem, result [, control])
```

The two names are compatibility entry points to the same translated solver.
Analytic derivatives are used when their procedure pointers are associated.

## KKT diagnostics

```fortran
call kkt_diagnose(problem, result, diagnostics [, tol])
```

The returned `kkt_diagnostics` contains stationarity, equality, inequality,
bound, dual-feasibility, and complementarity measures plus logical
`primal_feasible` and `first_order` flags.

## Starting points and multistart

```fortran
call startpars(problem, n_starts, starts [, seed, include_start, &
               feasibility_iter, status, message])
call csolnp_ms(problem, n_starts, result [, control, seed])
call gosolnp(problem, n_starts, result [, control, seed])
```

`starts` has shape `(n_starts, problem%n)`. A Halton sequence makes generation
portable and repeatable. `csolnp_ms` and `gosolnp` are sequential. Results are
returned in `multistart_result`, including every run and the selected best run.

## Standardization

```fortran
call solnp_standardize_problem(problem, standardized [, status, message])
```

This records equality shifts and expands two-sided inequalities into one-sided
standard-form metadata. The solver itself accepts either original or
standardized problem objects.

## Benchmark suite

```fortran
call solnp_problems_table(entries)
call solnp_problem_suite(suite, number, problem [, status, message])
```

The registry contains 68 Hock-Schittkowski table entries and nine `Other`
entries. The executable definitions are:

- HS01, HS05, HS06, HS11, HS14, HS24, HS38, HS64, HS65
- alkylation, box, entropy, garch, himmelblau5, powell,
  rosen_suzuki, wright4, and wright9

Use suite `"hs"` or `"Hock-Schittkowski"` with the HS number, or suite
`"Other"` with numbers 1 through 9 in the order shown by
`solnp_problems_table`.
