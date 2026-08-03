# API

Import the public interface with:

```fortran
use maxlik
```

All real calculations use `dp = kind(1.0d0)`. Callback status value zero means
success; a nonzero callback status is converted to `MAXLIK_EVALUATION_ERROR`.

## Problem definition

### `type(maxlik_problem)`

Important components:

- `npar`, `nobs`
- procedure pointers `objective`, `gradient`, `hessian`, and `scores`
- logical vector `active`
- vectors `lower`, `upper`
- optional linear-constraint arrays `eq_a`, `eq_b`, `ineq_a`, `ineq_b`

### `initialize_problem(problem, npar, objective, nobs)`

Initializes dimensions, attaches the scalar objective callback, activates every
parameter, and supplies unbounded defaults. `nobs` is optional but required when
an observation-score callback will be used.

Attach optional callbacks directly:

```fortran
problem%gradient => gradient_callback
problem%hessian  => hessian_callback
problem%scores   => score_callback
```

The callback interfaces are:

```fortran
subroutine objective(x, value, status)
subroutine gradient(x, gradient, status)
subroutine hessian(x, hessian, status)
subroutine scores(x, scores, status)
```

The score matrix has shape `nobs x npar`; rows are observation contributions
whose column sums form the full gradient.

### Parameter and constraint setup

- `set_fixed(problem, indices)` marks integer parameter positions fixed.
- `clear_fixed(problem)` reactivates all parameters.
- `set_bounds(problem, lower, upper, status)` sets box bounds.
- `set_equality_constraints(problem, a, b, status)` sets `a*x + b = 0`.
- `set_inequality_constraints(problem, a, b, status)` sets `a*x + b >= 0`.
- `clear_constraints(problem)` removes both linear-constraint sets.
- `constraint_violation(problem, x)` returns the maximum absolute equality
  residual or negative inequality residual.

## Control

### `type(maxlik_control)`

Defaults follow the upstream package where practical. Major fields include:

- `tol`, `reltol`, `gradtol`, `steptol`
- `qrtol`, `qac`
- `marquardt_lambda0`, `marquardt_lambda_step`,
  `marquardt_max_lambda`
- `nm_alpha`, `nm_beta`, `nm_gamma`
- `sann_temp`, `sann_tmax`, `random_seed`
- `learning_rate`, `batch_size`, `gradient_clip`
- `sga_momentum`, `adam_momentum1`, `adam_momentum2`
- `patience`, `patience_step`
- `iterlim`, `final_hessian`
- `store_values`, `store_parameters`
- `constraint_max_outer`, `constraint_tol`, `constraint_rho0`,
  `constraint_rho_factor`

## Maximization

### `max_lik(problem, start, result, method, control)`

`method` is optional and defaults to `nr`. Accepted names are:

- `nr`, `newton`, `newton-raphson`
- `bfgs`
- `bfgsr`, `bfgs-r`
- `bhhh`
- `cg`, `conjugate-gradient`
- `nm`, `nelder-mead`
- `sann`, `simulated-annealing`
- `sga`, `stochastic-gradient`
- `adam`

Direct solver entry points are also public: `solve_newton`, `solve_bfgs`,
`solve_bfgsr`, `solve_bhhh`, `solve_cg`, `solve_nelder_mead`, `solve_sann`,
`solve_sga`, and `solve_adam`.

BHHH, SGA, and Adam require `problem%scores` and a positive `problem%nobs`.
The deterministic methods can use analytic derivatives or numerical fallbacks.
Nelder-Mead and simulated annealing use only objective values.

## Results

### `type(maxlik_result)`

Important components:

- `estimate`, `maximum`
- `gradient`, `hessian`
- `covariance`, `std_error`, `condition_number`
- `gradient_obs`
- `active`
- `code`, `message`, `converged`, `method`
- `iterations`, `outer_iterations`
- `function_count`, `gradient_count`, `hessian_count`
- `constraint_violation`
- `stored_values`, `stored_parameters`

Convergence codes preserve the important upstream meanings:

- `MAXLIK_SUCCESS_GRADIENT = 1`
- `MAXLIK_SUCCESS_VALUE = 2`
- `MAXLIK_STEP_FAILURE = 3`
- `MAXLIK_ITERATION_LIMIT = 4`
- `MAXLIK_CONSTRAINT_FAILURE = 5`
- `MAXLIK_INVALID_START = 100`

Additional explicit input/evaluation/singularity codes are also exported.

## Numerical derivatives and utilities

- `numeric_gradient(problem, x, gradient, function_count, status, central)`
- `numeric_hessian(problem, x, hessian, function_count, gradient_count, status)`
- `compare_derivatives(problem, x, comparison, tolerance)`
- `numeric_jacobian(function, x, number_values, jacobian, status, step, active)`
- `pack_active_parameters` and `unpack_active_parameters`
- `progressive_condition_numbers(matrix, values, status, normalize)`

`type(derivative_comparison)` stores analytic and numerical arrays, errors,
maximum errors, and pass/fail flags.

## Inference

- `covariance_matrix(hessian, active, covariance, status)` returns the inverse
  observed information on active parameters.
- `robust_covariance_matrix(hessian, scores, active, covariance, status)`
  computes the score-sandwich covariance.
- `standard_errors(covariance, std_error)`
- `normal_confidence_intervals(estimate, std_error, z_value, lower, upper)`
- `maxlik_aic(maximum, number_parameters)`
- `condition_number(hessian, active, status)`

Fixed-parameter covariance rows and columns are zero.
