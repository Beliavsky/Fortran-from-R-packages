# API mapping

| R `tabuSearch` API | Fortran API | Status |
|---|---|---|
| `tabuSearch(...)` | `run_tabu_search(size, objective, result, control, initial_config)` | Implemented |
| `summary.tabu()` numerical calculations | `summarize_tabu(result, summary)` | Implemented |
| `plot.tabu()` | none | Intentionally omitted (plotting) |

## Search arguments

| R argument | Fortran equivalent |
|---|---|
| `size` | first argument of `run_tabu_search` |
| `iters` | `tabu_control%iters` |
| `objFunc` | Fortran procedure callback `tabu_objective` |
| `config` | optional `initial_config` |
| `neigh` | `tabu_control%neigh`; zero means `size` |
| `listSize` | `tabu_control%list_size` |
| `nRestarts` | `tabu_control%n_restarts` |
| `repeatAll` | `tabu_control%repeat_all` |
| `verbose` | `tabu_control%verbose` |

## Result mapping

| R result field | Fortran field/method |
|---|---|
| `configKeep` | `tabu_result%config_keep` |
| `eUtilityKeep` | `tabu_result%utility_keep` |
| `iters` | `tabu_result%iters` |
| `neigh` | `tabu_result%neigh` |
| `listSize` | `tabu_result%list_size` |
| `repeatAll` | `tabu_result%repeat_all` |
| maximum objective | `tabu_result%best_value()` |
| optimum configuration | `tabu_result%best_configuration()` |

`tabu_summary_result` additionally contains the optimum iterations, optimum
number of selected variables, all optimum occurrences, per-variable selection
counts, move frequencies, and the number of unique configurations visited.
