# API reference

All public numerical interfaces are available through:

```fortran
use risk_parity_portfolio_mod
```

Real calculations use `real(dp)`, where `dp = kind(1.0d0)`.

## High-level routine

### `risk_parity_portfolio`

```fortran
call risk_parity_portfolio(sigma, result, b, mu, lambda_mu, lambda_var, &
                           lower, upper, cmat, cvec, dmat, dvec, &
                           formulation, method_init, w0, theta0, gamma, &
                           zeta, tau, maxiter, ftol, wtol)
```

Only `sigma` and `result` are required.

- `sigma(n,n)`: symmetric covariance matrix with positive diagonal.
- `b(n)`: positive risk-budget targets; normalized internally.
- `mu(n)`: expected returns.
- `lambda_mu`: reward on expected return.
- `lambda_var`: penalty on portfolio variance.
- `lower(n)`, `upper(n)`: weight bounds; defaults are zero and one.
- `cmat(m,n)*w = cvec(m)`: additional equalities.
- `dmat(k,n)*w <= dvec(k)`: additional inequalities.
- `formulation`: one of the constants below. If absent for the vanilla
  long-only problem, a convex fast solver is used.
- `method_init`: initial vanilla solver constant.
- `w0(n)`: initial portfolio for SCA.
- `theta0`: initial auxiliary theta for theta formulations.
- `gamma`, `zeta`, `tau`: SCA learning-rate, decay, and regularization values.
- `maxiter`, `ftol`, `wtol`: iteration and convergence controls.

The routine always appends `sum(w)=1` and appends the supplied box bounds as
linear inequalities.

## Result type

`type(risk_parity_result)` contains:

- `weights(:)`
- `relative_risk_contribution(:)`
- `objective_history(:)`
- `risk_concentration`
- `mean_return`
- `variance`
- `theta`
- `iterations`
- `status`
- `converged`
- `feasible`
- `method`
- `formulation`

Status constants are `RPP_OK`, `RPP_INVALID_INPUT`, `RPP_INFEASIBLE`,
`RPP_LINEAR_SOLVE_FAILED`, and `RPP_MAX_ITER`.

## Initial-solver constants

- `INIT_CYCLICAL_SPINU`
- `INIT_CYCLICAL_RONCALLI`
- `INIT_CYCLICAL_CHOI`
- `INIT_NEWTON`

## Formulation constants

- `FORM_RC_DOUBLE_INDEX`
- `FORM_RC_OVER_B_DOUBLE_INDEX`
- `FORM_RC_OVER_VAR_VS_B`
- `FORM_RC_OVER_VAR`
- `FORM_RC_OVER_SD_VS_B_SD`
- `FORM_RC_VS_B_VAR`
- `FORM_RC_VS_THETA`
- `FORM_RC_OVER_B_VS_THETA`

The first, fourth, and seventh formulations target equal contributions rather
than arbitrary `b`. The other formulations use `b` explicitly.

## Direct solvers

```fortran
call risk_parity_ccd_spinu(sigma, b, w, info, iterations, tol, maxiter)
call risk_parity_ccd_roncalli(sigma, b, w, info, iterations, tol, maxiter)
call risk_parity_ccd_choi(sigma, b, w, info, iterations, tol, maxiter)
call risk_parity_newton(sigma, b, w, info, iterations, tol, maxiter)
call active_risk_parity_ccd(sigma, b, mu, tradeoff, risk_free, w, &
                            info, iterations, tol, maxiter)
```

## SCA routine

`risk_parity_sca` accepts complete constraint matrices, including the budget
constraint, and is useful when an application constructs its own constraint
system. Most users should call `risk_parity_portfolio`.

## Risk quantities and objectives

- `portfolio_variance(sigma,w)`
- `relative_risk_contributions(sigma,w)`
- `diagonal_risk_parity(sigma,b)`
- `objective_spinu(sigma,x,b)`
- `objective_roncalli(sigma,x,b)`
- `risk_vector(x,sigma,b,formulation)`
- `risk_jacobian(x,sigma,b,formulation)`
- `risk_objective(x,sigma,b,formulation)`
- `risk_gradient(x,sigma,b,formulation)`

## Projection and QP utilities

- `project_budget_box`
- `project_equality`
- `project_feasible`
- `project_linear_constraints`
- `is_feasible`
- `solve_equality_qp`
- `solve_qp_active_set`

The QP convention is:

```text
minimize  0.5*x^T*Q*x + q^T*x
subject to C*x = c
           D*x <= d
```
