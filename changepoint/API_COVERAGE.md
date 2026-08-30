# API coverage

This file records computational coverage relative to R package `changepoint`
2.3. The Fortran interface is intentionally typed and array-oriented rather
than an emulation of R's S4 classes.

| Upstream area | Fortran coverage |
| --- | --- |
| `PELT` and distribution-specific multiple-change wrappers | `cp_pelt` with the six `cp_cost_*` models; pruning is retained |
| `BINSEG` and distribution-specific wrappers | `cp_binseg` |
| `SEGNEIGH` and distribution-specific wrappers | `cp_segneigh`, a unified exact segment-neighborhood dynamic program |
| Normal mean | `cp_cost_mean_normal` |
| Normal variance with common/known mean | `cp_cost_var_normal`, optional `known_mean` |
| Normal mean and variance | `cp_cost_meanvar_normal` |
| Exponential | `cp_cost_exponential` |
| Gamma | `cp_cost_gamma`, optional/required `shape` as appropriate |
| Poisson | `cp_cost_poisson` |
| AMOC likelihood methods | `cp_amoc` |
| Mean CUSUM | `cp_amoc_cusum`, `cp_binseg_cusum`, `cp_segneigh_cusum` |
| Variance CSS | `cp_amoc_css`, `cp_binseg_css`, `cp_segneigh_css` |
| CROPS / range of penalties | `cp_crops` returning explicit penalty intervals and segmentations |
| `cpt.reg` Normal regression changes | `cp_regression_amoc`, `cp_regression_pelt`, `cp_regression_segment_fit` |
| `penalty_decision` | `cp_penalty_value` for BIC/SIC, AIC, Hannan-Quinn, MBIC, none, manual, and supported asymptotic penalties |
| `decision` | `cp_decision` |
| AMOC class-free asymptotic confidence transforms | `cp_amoc_asymptotic_value` for the upstream-supported scalar distribution cases |
| `fit.mean` | `cp_segment_means` |
| `fit.var` | `cp_segment_variances_mle` |
| `fit.scale` | `cp_segment_scales` |
| `fit.trend` | `cp_segment_trend_fits` |
| `fit.reg` | `cp_segment_regression_fits` |

## Intentional differences and omissions

- S4 classes, replacement/accessor methods, `show`, `summary`, `plot`, `acf`,
  `fitted`, `residuals`, `logLik`, `zoo`/`ts` time-coordinate presentation,
  formula parsing, and matrix-of-datasets convenience dispatch are R interface
  facilities and are not reproduced.
- Manual penalties are numeric inputs. Arbitrary R character expressions are
  not parsed or evaluated.
- The native segment-neighborhood implementation is one exact dynamic program
  over the same segment costs instead of duplicating the separate historical R
  functions for each distribution. It also accepts an explicit `minseglen`.
- Native Binary Segmentation and regression routines interpret `minseglen`
  literally. They do not intentionally reproduce historical C off-by-one
  quirks that were changed across upstream releases.
- Exponential input follows upstream behavior and permits zero observations.
  An all-zero segment is guarded with the smallest positive representable scale
  before taking a logarithm, avoiding a non-finite `log(0)` result.
- Poisson input is validated as finite, nonnegative integer-valued data.
- `fit.meanar` and `fit.trendar` are legacy class-specific parameter helpers
  retained in upstream `fit.R`, but current public changepoint constructors do
  not expose corresponding model families. They are not part of this native
  computational API.
- The R package's object-construction and warning behavior is represented by
  integer status values and typed result objects rather than reproduced textually.

## Result types

`changepoint_result`, `amoc_result`, `binseg_result`, `segneigh_result`, and
`crops_solution` hold native outputs. Status constants are `cp_ok`,
`cp_invalid_argument`, `cp_invalid_data`, and `cp_linalg_failure`.
