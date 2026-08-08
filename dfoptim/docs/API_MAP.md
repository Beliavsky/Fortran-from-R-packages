# API crosswalk

| R routine | Fortran routine | Notes |
|---|---|---|
| `hjk(par, fn, control, ...)` | `hjk(par, fn, control, user_data, monitor)` | Unbounded Hooke-Jeeves |
| `hjkb(par, fn, lower, upper, control, ...)` | `hjkb(par, fn, lower, upper, control, user_data, monitor)` | Elementwise box bounds |
| `nmk(par, fn, control, ...)` | `nmk(par, fn, control, user_data, monitor)` | Modified Nelder-Mead with restarts |
| `nmkb(par, fn, lower, upper, control, ...)` | `nmkb(par, fn, lower, upper, control, user_data, monitor)` | Uses the same hyperbolic/log transformations |
| `mads(par, fn, lower, upper, scale, control, ...)` | `mads(par, fn, lower, upper, scale, control, user_data, monitor)` | Lower-triangular randomized polling |

## Control mapping

### Hooke-Jeeves

| R control | Fortran component |
|---|---|
| `tol` | `hj_control_t%tol` |
| `maxfeval` | `hj_control_t%maxfeval` |
| `maximize` | `hj_control_t%maximize` |
| `target` | `hj_control_t%target` |
| `info` | `hj_control_t%trace` |
| not exposed | `hj_control_t%seed` |

### Nelder-Mead

| R control | Fortran component |
|---|---|
| `tol` | `nmk_control_t%tol` |
| `maxfeval` | `nmk_control_t%maxfeval` |
| `regsimp` | `nmk_control_t%regular_simplex` |
| `maximize` | `nmk_control_t%maximize` |
| `restarts.max` | `nmk_control_t%max_restarts` |
| `trace` | `nmk_control_t%trace` |

### MADS

| R control | Fortran component |
|---|---|
| `trace` | `mads_control_t%trace` |
| `tol` | `mads_control_t%tol` |
| `maxfeval` | `mads_control_t%maxfeval` |
| `maximize` | `mads_control_t%maximize` |
| `pollStyle="lite"` | `mads_poll_lite` |
| `pollStyle="full"` | `mads_poll_full` |
| `deltaInit` | `mads_control_t%delta_init` |
| `expand` | `mads_control_t%expand` |
| `lineSearch` | `mads_control_t%line_search` |
| `seed` | `mads_control_t%seed` |

## R infrastructure not translated

R lists and data frames are represented by derived types and allocatable
arrays. R's `...`, `match.fun`, console warnings, and global random-number state
are replaced by explicit Fortran callbacks, status codes, and local RNG state.
