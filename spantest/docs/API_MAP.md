# API map

| R procedure | Fortran procedure | Module | Notes |
|---|---|---|---|
| `span_bj` | `span_bj` | `spantest_classical` | `type(span_result)` |
| `span_f1` | `span_f1` | `spantest_classical` | `type(span_result)` |
| `span_f2` | `span_f2` | `spantest_classical` | `type(span_result)` |
| `span_grs` | `span_grs` | `spantest_classical` | ML residual covariance retained |
| `span_hk` | `span_hk` | `spantest_classical` | Joint frontier test |
| `span_km` | `span_km` | `spantest_classical` | Difference-return regression |
| `span_py` | `span_py` | `spantest_classical` | Normal-reference statistic |
| `span_gl_a` | `span_gl_a` | `spantest_gl` | Seed and simulation count are optional arguments |
| `span_gl_ad` | `span_gl_ad` | `spantest_gl` | Joint alpha/delta restrictions |
| `span_as` | `span_as` | `spantest_as` | Named p-values in `type(as_result)` |
| `span_simulate` | `span_simulate` | `spantest_simulation` | Typed matrices in `type(simulation_result)` |
| `garch_filter` | `garch_filter` | `spantest_simulation` | Public computational kernel |
| `f_cauchypv` | `cauchy_pvalue` | `spantest_as` | Cauchy p-value merger |
| `f_rsstd` | `standardized_skew_t` | `spantest_simulation` | Fernandez-Steel skew-t generator |
| `gl_sim_stats` | internal loop in `span_gl_*` | `spantest_gl` | Streaming implementation |

R lists are represented by derived types. R's `NA` results are represented by
IEEE quiet NaNs together with a nonzero status code and diagnostic message.
