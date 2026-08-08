# API map

## Callback classes

| RcppNumerical C++ | Fortran |
|---|---|
| `Func::operator()(double)` | `scalar_function_interface` |
| `MFunc::operator()(Constvec)` | `multivariate_function_interface` |
| `MFuncGrad::f_grad()` | `objective_gradient_interface` |
| C++ functor state | optional polymorphic `user_data` |

## Integration

| RcppNumerical API | Fortran API | Notes |
|---|---|---|
| scalar `integrate()` | `integrate_1d()` | Same default tolerances and rule 41 |
| `GaussKronrod15` ... `GaussKronrod201` | `gk15` ... `gk201` | All 12 rules retained |
| scalar output references | `integration_result_t` | Value, error, code, evaluations, intervals |
| multivariate `integrate()` | `integrate_nd()` | Default Cuhre rules 13/11/9 |
| Cuba failure flag | `cubature_result_t%error_code` | 0 success, 1 evaluation limit, negative input error |

`integrate_nd()` supports finite, semi-infinite, and doubly infinite bounds.
Dimensions 2 through 20 are supported. The upper bound must exceed the lower
bound in every dimension.

## Optimization

| RcppNumerical API | Fortran API |
|---|---|
| `optim_lbfgs()` | `optim_lbfgs()` |
| `optim_lbfgsb()` | `optim_lbfgsb()` |
| integer status and `fx_opt` | `optimization_result_t` |

The wrapper-level status remains 0 for convergence and -1 for failure. The
underlying optimizer status is retained in `native_status`.

## Logistic regression

| R API | Fortran API |
|---|---|
| `fastLR(x, y, start, eps_f, eps_g, maxit)` | `fast_lr(x, y, fit, start, eps_f, eps_g, maxit)` |
| R list | `logistic_fit_t` |

The returned type contains coefficients, fitted probabilities, linear
predictors, maximized log likelihood, iteration/evaluation counts, and a
convergence flag.
