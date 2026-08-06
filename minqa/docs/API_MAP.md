# API map

| R API | Modern Fortran API | Coverage |
|---|---|---|
| `bobyqa(par, fn, lower, upper, control, ...)` | `call bobyqa(objective, x, result, lower, upper, control)` | Complete computational translation |
| `newuoa(par, fn, control, ...)` | `call newuoa(objective, x, result, control)` | Complete computational translation |
| `uobyqa(par, fn, control, ...)` | `call uobyqa(objective, x, result, control)` | Complete computational translation |
| `commonArgs` | `minqa_control_t` plus native validation | Replaced by typed control handling |
| `print.minqa` | User formats `minqa_result_t` directly | R-only presentation layer omitted |
| R closure and `...` forwarding | Fortran procedure callback | Replaced by native procedure interface |
| Rcpp `.Call` wrappers | Direct module calls | Removed |
| R evaluation environment/counter | `result%evaluations` | Replaced natively |

## Control mapping

| R control entry | Fortran component |
|---|---|
| `npt` | `control%npt` |
| `rhobeg` | `control%rhobeg` |
| `rhoend` | `control%rhoend` |
| `iprint` | `control%iprint` |
| `maxfun` | `control%maxfun` |
| `obstop = FALSE` | `control%adjust_start = .true.` |
| `force.start` | Not needed as a separate option; BOBYQA core adjustment is retained |

The Fortran result contains the optimized vector, objective value,
evaluation count, mapped status, raw status, and message.
