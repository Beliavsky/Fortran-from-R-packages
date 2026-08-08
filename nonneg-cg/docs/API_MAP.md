# API map

| R/C interface | Fortran interface | Notes |
|---|---|---|
| `minimize.nonneg.cg()` | `minimize_nonneg_cg()` | Native typed callback interface |
| `evaluate_function` | `objective_callback` | Scalar objective callback |
| `evaluate_gradient` | `gradient_callback` | Analytic gradient callback |
| `x0` | `x` | Must be feasible; updated in place |
| `tol` | `nonneg_cg_control_t%tol` | Same stopping quantity, `abs(dot(g,d))` |
| `maxnfeval` | `%maxnfeval` | Original counter semantics preserved |
| `maxiter` | `%maxiter` | Zero means unlimited |
| `decr_lnsrch` | `%decr_lnsrch` | Same geometric backtracking factor |
| `lnsrch_const` | `%lnsrch_const` | Same sufficient-decrease constant |
| `max_ls` | `%max_ls` | Same line-search trial limit |
| `extra_nonneg_tol` | `%extra_nonneg_tol` | Explicitly clamps nonpositive values to zero |
| R result `x` | `result%x` | Final parameters |
| R result `fun` | `result%fun` | Final objective |
| R result `niter` | `result%niter` | Same iteration convention |
| R result `nfeval` | `result%nfeval` | Same legacy counter convention |
| R result `term` | `result%status`, `result%message` | Integer constants and text |

The Rcpp wrapper, dynamic registration, R lists, and `do.call` machinery are
interface code and are not translated.
