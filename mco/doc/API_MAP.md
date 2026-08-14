# API map

| R `mco` interface | Fortran interface | Coverage |
|---|---|---|
| `nsga2()` | `nsga2_optimize()` | Native callback-based NSGA-II |
| `paretoFilter.matrix()` | `pareto_filter()` / `pareto_mask()` | Direct numerical equivalent |
| `paretoSet()` | `result%par(:, result%pareto_optimal)` | Native result component |
| `paretoFront()` | `result%value(:, result%pareto_optimal)` | Native result component |
| `normalizeFront()` | `normalize_front()` | Direct equivalent |
| `generationalDistance()` | `generational_distance()` | Direct equivalent |
| `generalizedSpread()` | `generalized_spread()` | Standard absolute-deviation form |
| `dominatedHypervolume()` | `dominated_hypervolume()` | Exact recursive slicing |
| `epsilonIndicator()` | `epsilon_indicator()` | Upstream additive indicator |
| exported test functions | same lower-case procedure names | Direct equations |
| plot methods | omitted | Noncomputational |
| S3 result classes | `nsga2_result` derived type | Native Fortran replacement |
| vectorized R callbacks | scalar Fortran callbacks | Omitted adapter layer |
| vector-valued generation snapshots | final population result | Omitted |

All fronts are stored as `objective by point` matrices in Fortran.
