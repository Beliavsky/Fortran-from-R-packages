# API map

| R API / subsystem | Fortran API | Coverage |
|---|---|---|
| `mize()` | `mize_minimize` | Native one-shot optimizer |
| `make_mize()` | `mize_control_t` plus `mize_state_t` | Typed replacement |
| `mize_init()` | `mize_init` | Implemented |
| `mize_step()` | `mize_step` | Implemented |
| `mize_step_summary()` | Public `mize_state_t` fields and `mize_state_result` | Typed replacement |
| `check_mize_convergence()` | `check_mize_convergence` | Implemented |
| `check_mize_gradient()` | `check_mize_gradient`, `gradient_check_t` | Implemented |
| SD | `control%method = 'SD'` | Implemented |
| BFGS | `control%method = 'BFGS'` | Implemented |
| SR1 | `control%method = 'SR1'` | Implemented |
| L-BFGS | `control%method = 'L-BFGS'` | Implemented |
| nonlinear CG | `control%method = 'CG'` | All 11 update formulas implemented |
| Newton | `control%method = 'Newton'` | Exact callback or finite-difference Hessian |
| partial Hessian | `control%method = 'pHess'` | Implemented |
| truncated Newton | `control%method = 'TN'` | Hessian-vector or finite-difference Newton-CG |
| NAG | `control%method = 'NAG'` | Implemented with native look-ahead gradient |
| momentum | `control%method = 'MOM'` | Classical/Nesterov styles and schedules |
| delta-bar-delta | `control%method = 'DBD'` | Implemented |
| constant step | `line_search = 'constant'` | Implemented |
| backtracking / Armijo | `line_search = 'backtracking'` | Implemented |
| bold driver | `line_search = 'bold driver'` | Implemented |
| More-Thuente, Rasmussen, Schmidt, Hager-Zhang | Corresponding `line_search` strings | Routed through the shared safeguarded Wolfe engine; condition defaults remain configurable |
| R lifecycle hooks and dependency injection | Monitor and momentum callbacks | Numerical callback replacement, not a general aspect system |
| R progress data frame | Allocatable progress arrays in `mize_result_t` | Implemented |

The Fortran API uses a combined objective-gradient callback because it avoids
R-style list dispatch and permits a single evaluation to share intermediate
work. Function and gradient evaluation counts therefore advance together for
that callback.
