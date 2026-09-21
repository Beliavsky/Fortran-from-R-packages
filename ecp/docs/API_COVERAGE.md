# API coverage

The coverage basis is the nine functions explicitly exported by upstream `NAMESPACE`. All nine are computational and all nine have meaningful Fortran mappings.

| R export | Fortran API | Main compatibility notes |
| --- | --- | --- |
| `e.agglo` | `e_agglo` | Cyclic energy-distance merging is retained. An arbitrary R penalty function cannot be represented directly; callers may provide stage penalty values instead. |
| `e.divisive` | `e_divisive` | Energy splits and within-segment permutation tests are retained. Fixed-`k` mode is deterministic; permutation mode uses `random_number`, not R RNG state. |
| `getBetween` | `get_between` | Direct numerical translation. |
| `getWithin` | `get_within` | Direct numerical translation. |
| `e.cp3o` | `e_cp3o` | Energy objective, DP paths, candidate pruning, and cached/prefix-sum statistics are translated; `prune=.false.` retains exhaustive verification. |
| `e.cp3o_delta` | `e_cp3o_delta` | Windowed-energy objective, pruning, and DP paths are translated; prefix sums evaluate window/tail terms and `prune=.false.` retains exhaustive verification. |
| `ks.cp3o` | `ks_cp3o` | Scalar KS behavior, pruning, cached candidate statistics, and DP paths are translated. Multivariate input uses maximum marginal KS rather than upstream matrix linear indexing. |
| `ks.cp3o_delta` | `ks_cp3o_delta` | Windowed KS pruning and DP paths are translated; maximum marginal rule for multivariate data; exhaustive verification is available with `prune=.false.`. |
| `kcpa` | `kcpa` | Kernel costs, DP and model-selection penalty are translated; O(n^2) prefix-sum segment costs replace direct block summation, and optional `bandwidth_rows` can reproduce a chosen upstream subsample. |

Internal upstream helpers such as `splitPoint`, `sig.test`, `process.data`, `srcGetV`, and `srcKcpa` are not separate coverage entries because the requested denominator counts distinct exported computational R functions. Their numerical operations are implemented as internal Fortran procedures where needed.
