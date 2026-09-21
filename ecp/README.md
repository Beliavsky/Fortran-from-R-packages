# ecp - modern Fortran translation

This directory translates the computational surface of the R package **ecp 3.1.6** to modern free-form Fortran with FPM. The focus is the package's nonparametric multiple change-point algorithms: energy statistics, hierarchical splitting/merging, CP3O, Kolmogorov-Smirnov CP3O, and kernel change-point analysis.

The maintained implementation is pure Fortran. It does not vendor Rcpp, BLAS, LAPACK, or another translated R package. A single package-local `dp = real64` kind is defined in `ecp_kinds` and used throughout.

## Build

```text
fpm build
fpm test
fpm run --example change_points
```

The public Fortran facade is `ecp_api`. Important result types are `cp_result`, `agglo_result`, and `divisive_result`.

## Translation coverage

Package status: **substantial**.

**9 of 9 (100.0%)** exported computational R functions are mapped. This percentage measures mapped computational functions under the requested coverage convention; it does **not** mean complete R compatibility. In particular, R list/S3 presentation, timing/verbose output, arbitrary penalty closures, and exact R RNG streams are not part of the typed Fortran interface. CP3O candidate pruning is implemented and can be disabled with `prune=.false.` for exhaustive verification.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `e.agglo` | `ecp_agglo` | `e_agglo` | substantial |
| `e.divisive` | `ecp_divisive` | `e_divisive` | substantial |
| `getBetween` | `ecp_energy` | `get_between` | complete |
| `getWithin` | `ecp_energy` | `get_within` | complete |
| `e.cp3o` | `ecp_cp3o` | `e_cp3o` | substantial |
| `e.cp3o_delta` | `ecp_cp3o` | `e_cp3o_delta` | substantial |
| `ks.cp3o` | `ecp_cp3o` | `ks_cp3o` | substantial |
| `ks.cp3o_delta` | `ecp_cp3o` | `ks_cp3o_delta` | substantial |
| `kcpa` | `ecp_kernel` | `kcpa` | substantial |

See `docs/API_COVERAGE.md` for interface differences and `docs/PROVENANCE.md` for upstream provenance.

## Algorithm notes

`e_cp3o` and the three related CP3O procedures use the upstream adjacent-segment objective, dynamic-programming state, and candidate-pruning inequality. An optional `prune=.false.` argument retains the exhaustive search as a built-in reference path. `cp_result` reports `statistic_evaluations` and `candidates_pruned` for diagnostics. For energy CP3O, a two-dimensional prefix sum of the pairwise distance matrix makes each full energy statistic O(1) after the O(n^2) distance preprocessing; the delta-energy path also uses an adjacent-distance prefix. Pruning statistics are cached and reused by the next DP level.

The upstream KS C++ code indexes an R matrix through one-dimensional `NumericMatrix::operator[]` expressions in several places. The Fortran translation matches the scalar case directly and defines multivariate behavior as the maximum marginal KS statistic. The same pruning framework is used for both KS variants, with retained statistics cached between levels.

For `kcpa`, the Gaussian-kernel cost and dynamic program are retained. Kernel block-prefix sums reduce construction of all segment costs from O(n^4) direct block summation to O(n^2). The default bandwidth calculation remains deterministic and uses all observations. Callers that want to reproduce the upstream `n > 250` sampling choice can pass the chosen 1-based rows through `bandwidth_rows`; exact R RNG streams are intentionally not reproduced.

## License

The upstream package declares `GPL (>= 2)`. This translation is distributed under GPL-2.0-or-later. See `LICENSE`, `NOTICE.md`, and `docs/PROVENANCE.md`.
