# API map

| Upstream R code | Modern Fortran API | Notes |
|---|---|---|
| `solnl()` | `solnl()` | Main SQP solver; vector variables replace R matrices. |
| `objfun` callback | `objective_function` | Scalar procedure callback. |
| `confun` callback | `constraint_function` | Returns allocatable equality and inequality vectors. |
| finite-difference objective gradient | internal derivative engine | Forward differences by default; central optional. |
| finite-difference constraint Jacobian | internal derivative engine | Analytic callback also supported. |
| BFGS Hessian update | `damped_bfgs_update` | Powell-damped positive-definite update. |
| `quadprog::solve.QP()` | `quadprog::solve_qp()` | Supplied modern Fortran Goldfarb-Idnani solver. |
| internal `solqp()` fallback | elastic feasibility QP | A single nonnegative relaxation variable makes inconsistent linearizations solvable. |
| R list of multipliers | `type(nlc_multipliers)` | Linear/nonlinear equalities, bounds, and inequalities are separated. |
| R result list | `type(nlc_result)` | Adds status, message, constraint violation, KKT error, and QP counters. |
| `boundflag()` | bound-aware finite differences | Chooses forward/backward perturbations within finite bounds. |

The R-only matrix-shape restoration, printed warnings, list construction, and
exception handling are replaced by vector inputs and explicit status values.
