# Translation coverage

The R namespace exports many names because each numerical operation has formula,
model-class, print, plot, summary, and accessor methods. The Fortran port groups
these around explicit numerical entry points.

| Upstream area | Fortran coverage |
|---|---|
| `segmented.lm`, `seg.lm.fit`, `seg.num.fit` | `fit_segmented_lm`, `segmented_lm`, `segmented_numeric` |
| `segmented.glm`, `seg.glm.fit` | `fit_segmented_glm`, `segmented_glm` for Gaussian/binomial/Poisson |
| `stepmented.lm`, `step.lm.fit`, numeric variants | `fit_stepmented_lm`, `stepmented_lm`, `stepmented_numeric` |
| `stepmented.glm`, `step.glm.fit` | `fit_stepmented_glm`, `stepmented_glm` |
| `segmented.lme` | `fit_segmented_lme`, `segmented_lme`, backed by translated `nlme` |
| `segreg`, `stepreg` | Array-based `segreg` and `stepreg` dispatchers |
| `seg`, model-matrix helpers | `hinge_matrix`, `step_matrix`, `segmented_model_matrix` |
| `predict.segmented`, `predict.stepmented` | `predict_segmented` |
| `slope`, `intercept`, `broken.line`, `aapc` | `segment_slopes`, `segment_intercepts`, `broken_line_values`, `aapc` |
| `confint.segmented` | Delta/normal `breakpoint_confint` |
| `davies.test`, `pscore.test` | Grid Davies-style and fixed-breakpoint score tests |
| `pwr.seg` | Fixed-breakpoint score-test power approximation |
| `selgmented` | `select_breakpoints_bic` |
| summaries, coefficients, fitted values, covariance | Fields of `segmented_result` and `segmented_lme_result` |
| plotting, printing, drawing history | Omitted |
| bootstrap/reboot fitting | Not reproduced exactly |
| constrained segmented regression | Not reproduced exactly |
| segmented ARIMA | Not reproduced |
| stepmented time-series formula wrapper | Numeric step fitting is available; R time-series metadata is omitted |

All original R files are retained under `original/R`, including omitted methods.
