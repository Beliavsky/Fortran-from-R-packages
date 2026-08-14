# API mapping

## Exported computational API

| R API | Fortran API | Status |
|---|---|---|
| `kofnGA(...)` | `kofn_ga(n,k,objective,result,control,initpop)` | Implemented |
| `mutfrac` conversion | `mutation_probability(k,mutfrac)` | Implemented |
| `summary.GAsearch` numerical summary | `summarize_result(result)` | Implemented |

## Result fields

| R result | Fortran result |
|---|---|
| `bestsol` | `result%bestsol` |
| `bestobj` | `result%bestobj` |
| `pop` | `result%pop` |
| `obj` | `result%obj` |
| `old$best` | `result%best_history` |
| `old$obj` | `result%obj_history` |
| `old$avg` | `result%avg_history` |

## Intentionally omitted R infrastructure

- `plot.GAsearch`: plotting only.
- `print.GAsearch` and `print.summary.GAsearch`: display formatting only.
- `cluster`, `sharedmemory`, and helpers in `Rdsm_functions.R`: R process and
  shared-memory dispatch infrastructure. Objective evaluation remains a normal
  procedure callback in Fortran.
- R S3 class machinery.

The `bigmemory` package is therefore not a Fortran dependency.
