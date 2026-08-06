# API map

| Upstream R API | Fortran API | Notes |
|---|---|---|
| `tv` | `tv_transition` | Scalar speed, vector locations, vector transition variable |
| `tvObj` | `tv_objective_value` | Total and optional observationwise negative log likelihood |
| `garchObj` | `garchx_objective_value`, `garchx_filter` | Reused from supplied dependency |
| `tvgarchObj` | `tvgarch_objective_value` | Joint TV/GARCH objective |
| `tvgarchSim` | `tvgarch_simulate` | Typed simulation result |
| `tvgarch` | `fit_tvgarch` | Maximization by parts; optional joint refinement |
| `fitted.tvgarch` | `fitted_tvgarch` | `tv`, `garch`, or total variance |
| `predict.tvgarch` | `tvgarch_forecast` | Monte Carlo short-term forecast times future TV component |
| `quantile.tvgarch` | `tvgarch_quantile_path` | Empirical standardized-residual quantiles |
| `residuals.tvgarch` | `fit%residuals` | Direct typed component |
| `coef.tvgarch` | `fit%par_g`, `fit%hfit%par` | Direct typed components |
| `vcov.tvgarch` | `fit%vcov_g`, `fit%hfit%vcov` | Robust block estimates |
| `logLik.tvgarch` | `fit%loglik` | Gaussian log likelihood |
| `dccObj` | `dcc_objective`, `dcc_filter` | Objective, paths, observationwise terms |
| `mtvgarchSim` | `mtvgarch_simulate` | CCC or DCC simulation |
| `mtvgarch` | `fit_mtvgarch` | Equation-by-equation margins, CCC/DCC |
| spillover branch of `mtvgarch` | `fit_mtvgarch_spillover` | Iterated lagged cross-variance regressors |
| `predict.mtvgarch` | `mtvgarch_forecast` | Marginal variance forecasts |
| `combinations`, `combos` | `combinations_binary` | Nonempty binary subsets; ragged R form omitted |
| `tvgarchTest` | `tvgarch_test` | Robust and nonrobust TR2/LM tables and selected order |

Plot, print, summary, `toLatex`, date/index, and other presentation methods have
no Fortran computational equivalent and are omitted.
