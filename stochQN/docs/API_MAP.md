# API map

## Core optimizer classes

| R/C API | Fortran API | Notes |
|---|---|---|
| `oLBFGS_free`, `initialize_oLBFGS` | `type(olbfgs_t)` | Initialize, advance, reset, query iteration and memory use. |
| `SQN_free`, `initialize_SQN` | `type(sqn_t)` | Supports Hessian-vector and gradient-difference correction pairs. |
| `adaQN_free`, `initialize_adaQN` | `type(adaqn_t)` | Supports empirical-Fisher or gradient-difference correction pairs, AdaGrad/RMSProp scaling, and validation rollback. |
| `run_oLBFGS_free`, `run_oLBFGS` | `olbfgs_t%advance` | Reverse communication through `stochqn_request_t`. |
| `run_SQN_free`, `run_SQN` | `sqn_t%advance` | Request can contain both evaluation point and Hessian-vector direction. |
| `run_adaQN_free`, `run_adaQN` | `adaqn_t%advance` | Function value is supplied only for function requests. |
| `update_gradient` | gradient argument to `%advance` | Explicit typed array rather than mutable R list state. |
| `update_hess_vec` | Hessian-vector argument to `sqn_t%advance` | Explicit typed array. |
| `update_fun` | function-value argument to `adaqn_t%advance` | Explicit scalar. |
| `get_curr_x` | caller-owned `x` array | Optimizers update the array in place. |
| `get_iteration_number` | `%get_iteration()` | Available on all optimizer types. |

## Request and status values

The original integer codes are retained:

- tasks 100 through 105
- iteration information 200 through 203
- update statuses `-1000`, `0`, and `1`

The helper functions `stochqn_task_name`, `stochqn_info_name`, and
`stochqn_status_name` return descriptive strings.

## Guided wrappers

| R API | Fortran API |
|---|---|
| `oLBFGS` | `optimize_olbfgs` |
| `SQN` | `optimize_sqn` |
| `adaQN` | `optimize_adaqn` |
| `partial_fit` | reverse-communication `%advance`, or a complete callback runner |
| `predict.stochQN_guided` | application-defined prediction routine |

The callback gradient receives the task code, so it can select a new batch, the
same batch, or a larger accumulated batch.

## Logistic regression

| R API | Fortran API |
|---|---|
| `logistic_loss` | `logistic_loss` |
| `logistic_grad` | `logistic_gradient` |
| `logistic_Hess_vec` | `logistic_hessian_vector` |
| `logistic_pred` | `logistic_predict_probability` |
| `stochastic.logistic.regression` | `type(stochastic_logistic_t)` |
| `partial_fit_logistic` | `stochastic_logistic_t%partial_fit` |
| `coef.stoch_logistic` | `%get_coefficients()` |
| `predict.stoch_logistic` | `%predict_probability()` and `%predict_class()` |

The Fortran logistic model accepts numeric matrices and binary numeric targets.
R formula parsing, factors, model matrices, and class labels are intentionally
omitted.
