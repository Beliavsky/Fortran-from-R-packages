# API

Public module: `optimflex`.

## Solvers

- `bfgs(start, objective, result [, gradient, hessian, control])`
- `l_bfgs_b(start, objective, result [, lower, upper, gradient, hessian, control])`
- `newton_raphson(start, objective, result [, gradient, hessian, control])`
- `modified_newton(start, objective, result [, gradient, hessian, control])`
- `gauss_newton(start, objective, result [, residual, gradient, hessian, jacobian, control])`
- `levenberg_marquardt(start, objective, result [, lower, upper, residual, gradient, hessian, gn_hessian, jacobian, control])`
- `dogleg(start, objective, result [, lower, upper, residual, gradient, hessian, jacobian, control])`
- `double_dogleg(start, objective, result [, lower, upper, residual, gradient, hessian, gn_hessian, jacobian, control])`

Use keyword arguments when supplying optional callbacks.

## Differentiation

- `fast_grad`
- `fast_hess`
- `fast_jac`
- `get_eps`

Methods are `diff_forward`, `diff_central`, and `diff_richardson`.

## Types

`optim_control` contains the union of the upstream control-list fields.
Method-specific constructors are:

- `bfgs_default_control()`
- `lbfgsb_default_control()`
- `newton_default_control()`
- `modified_newton_default_control()`
- `gauss_newton_default_control()`
- `lm_default_control()`
- `dogleg_default_control()`
- `double_dogleg_default_control()`

`optim_result` contains `par`, `objective`, `converged`, `status`, `iter`,
`cpu_time`, `elapsed_time`, `max_grad`, `hess_is_pd`, `hessian`,
`approx_hessian`, `approx_hinv`, `pred_dec`, and `pred_dec_avg` where relevant.
