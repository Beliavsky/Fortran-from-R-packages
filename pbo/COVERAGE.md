# Computational coverage

Source package: `pbo` 1.3.5.

## Translated

| R functionality | Fortran functionality | Status |
|---|---|---|
| `pbo` CSCV partitioning | `compute_pbo` | Complete |
| `utils::combn` use | `generate_combinations` | Complete |
| User performance function | `performance_function` callback | Complete |
| IS/OOS strategy scoring | `performance_is`, `performance_oos` | Complete |
| IS and OOS maximizers | `selected_is`, `selected_oos` | Complete |
| Average OOS rank | `oos_rank` | Complete |
| Normalized rank and logit | `omega_bar`, `lambda` | Complete |
| PBO estimate | `phi` | Complete |
| Selected IS/OOS pairs | `selected_pairs` | Complete |
| Legacy `lm(rn_pairs)` summary | `slope`, `intercept`, `adjusted_r2` | Compatible |
| Intended OOS degradation fit | `degradation_*` | Added corrected result |
| Probability below threshold | `below_threshold` | Complete |
| Selection dotplot data | `selection_frequencies` | Complete |
| CSCV/pairs/ranks plot data | fields of `pbo_result` | Complete |
| Dominance plot CDF data | `dominance_curve` | Complete |
| Upstream `SD2` difference | `sd2_difference` | Complete |
| Integrated CDF difference | `integrated_difference` | Added |
| Sharpe example metric | `sharpe_ratio` | Complete |
| Omega example metric | `omega_ratio` | Complete |

## Not compiled

- Lattice, latticeExtra, and grid graphics.
- S3 print, summary, and plot dispatch.
- `foreach` and `doParallel` adapters.
- HTML/R Markdown vignette rendering.
- R data-frame, expression, and formula infrastructure.

Parallel execution can be added by callers around independent CSCV cases or by
a future OpenMP implementation. The present library intentionally has no
parallel-runtime dependency.
