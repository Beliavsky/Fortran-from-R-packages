# API map

This map distinguishes implemented numerical procedures from excluded R
infrastructure.

| Original R function | Fortran equivalent | Status |
|---|---|---|
| `garchxSim` | `garchx_simulate` | Implemented and tested |
| `garchxRecursion` | `garchx_filter`, `garchx_filter_derivatives` | Implemented and tested |
| `garchxObjective` | `garchx_objective_value` | Both objective modes implemented and tested |
| `garchx` | `fit_garchx` | Numerical fitting implemented and tested |
| `garchxAvar` | `garchx_asymptotic_covariance` | Original ordinary mode implemented and tested |
| `coef.garchx` | `garchx_fit%par` | Direct result field |
| `fitted.garchx` | `garchx_fit%sigma2` | Direct result field |
| `residuals.garchx` | `garchx_fit%residuals` | Direct result field |
| `logLik.garchx` | `garchx_fit%loglik`, `garchx_loglikelihood` | Implemented and tested |
| `nobs.garchx` | Array size after `garchx_max_lag` removal | Direct numerical operation |
| `vcov.garchx` | `garchx_covariance`, `garchx_fit%vcov` | Ordinary, robust, and HAC tested |
| `confint.garchx` | `confidence_intervals` | Implemented and tested |
| `predict.garchx` | `garchx_forecast` | Residual-bootstrap forecasts tested |
| `quantile.garchx` | `garchx_quantile_path` | Implemented and tested |
| `refit.garchx` | `refit_garchx` | Fixed and re-estimated modes tested |
| `glag` | Generic `glag` | Vector and matrix forms tested |
| `gdiff` | Generic `gdiff` | Vector and matrix forms tested |
| `rmnorm` | `rmnorm` | Implemented and moment-tested |
| `ttest0` | `boundary_t_tests` | Implemented and tested |
| `waldtest0` | `boundary_wald_test` | Implemented and simulation-tested |
| `print.garchx` | None | R output infrastructure excluded |
| `toLatex.garchx` | None | Formatting infrastructure excluded |
| S3 generic/method registration | None | R class infrastructure excluded |
| `zoo` indexing and trimming | None | Caller supplies clean numeric arrays |

The original `garchxAvar` source explicitly stops for robust and HAC covariance,
so those modes are not claimed for `garchx_asymptotic_covariance`.
