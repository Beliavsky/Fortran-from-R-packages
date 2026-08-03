# API map

## Primary routine

| Upstream R/C++ | Modern Fortran | Coverage |
|---|---|---|
| `GenSA(par, fn, lower, upper, control, ...)` | `gensa_minimize(objective, lower, upper, result, control, initial, constraint)` | Direct computational counterpart |
| Returned `value` | `result%value` | Direct |
| Returned `par` | `result%par` | Direct |
| Returned `counts` | `result%counts` | Direct |
| Returned `trace.mat` | `result%trace` | Same four quantities, recorded once per outer iteration |

## Control mapping

| R control | Fortran component | Notes |
|---|---|---|
| `maxit` | `maxit` | Maximum outer GSA iterations |
| `threshold.stop` | `has_threshold`, `threshold_stop` | The logical distinguishes an absent threshold from a numeric value |
| `nb.stop.improvement` | `no_improvement_stop` | Negative disables the test |
| `smooth` | `smooth` | Selects projected BFGS or pattern search |
| `max.call` | `max_call` | Counts all objective evaluations, including local search and numerical gradients |
| `max.time` | `max_time` | CPU seconds, matching the upstream `clock()` behavior |
| `temperature` | `temperature` | Initial visiting temperature |
| `visiting.param` | `visiting_param` | Must satisfy `1 < qv < 3` |
| `acceptance.param` | `acceptance_param` | Must be less than 1 |
| `verbose` | `verbose` | Periodic console reporting |
| `simple.function` | `simple_function` | Shortens the stagnation interval before local search |
| `trace.mat` | `trace` | Enables allocated trace arrays |
| `seed` | `seed` | Reinitializes the translated `ran2` generator |
| hidden `markov.length` | `markov_length` | Zero selects `2*dimension`; a supplied value must be a multiple of dimension |
| hidden `tem.restart` | `temp_restart` | Re-annealing threshold |
| hidden `high.dim` | `local_search` | Controls hybrid local searches directly |
| `REPORT` | `report` | Reporting interval when verbose |

Additional Fortran controls are `local_maxit`, `local_tolerance`, and
`max_constraint_attempts`.

## Internal numerical mapping

| Upstream routine | Fortran counterpart |
|---|---|
| `Engine::visita` | `gensa_rng::gensa_visit` |
| `Utils::ran2` | type-bound `ran2_state%uniform` |
| `Utils::yyGas` | type-bound `ran2_state%normal` |
| smooth `L-BFGS-B` local phase | `gensa_local::projected_bfgs` |
| hard `constrOptim` local phase | `gensa_local::bounded_pattern_search` |
| `Tracer` | `gensa_trace` arrays |

## Omitted R infrastructure

The `.Call` bridge, environments, names, list construction, `write.table`, R
warnings, and `...` argument dispatch are not reproduced. They do not alter
the numerical optimization model.
