# Source coverage

The table records how each upstream R source area is represented. Presentation
and R runtime adapters are intentionally excluded from the computational port.

| Upstream source | Fortran coverage |
|---|---|
| `R/specification.R` | `garch_spec`, `validate_specification`, model/distribution identifiers, initialization controls, variance targeting, and variance-regressor settings. |
| `R/constraints.R` | Model-specific bounds, reduced IGARCH/EWMA parameterizations, and stationarity/objective restrictions in `tsgarch_fit.f90`. |
| `R/initialization.R` | `initialize_parameters`, backcast/sample/unconditional variance initialization, and distribution defaults. The dormant unexported real-GARCH initializer is not included because `valid_garch_models()` excludes real-GARCH. |
| `R/likelihood.R` | `filter_garch` and `garch_loglikelihood`; distribution log densities come from vendored `tsdistributions`. |
| `R/filter.R` | Conditional variance, residual, component, and log-likelihood recursion in `tsgarch_model.f90`. |
| `R/estimate.R` | `estimate_garch`, parameter packing, bounded optimization, fit statistics, and optional inference. |
| `R/solvers.R` | Replaced by GPL-2-compatible bounded Nelder-Mead in `tsd_optimize.f90`; nloptr-specific option lists are not reproduced. |
| `R/distribution.R` | Ten standardized innovation laws and distribution moments supplied by vendored `tsdistributions` modules. |
| `R/simulate.R` and `src/simulation.cpp` | `simulate_garch` and `simulate_conditional`, including all valid models and distributions. Parametric simulation is implemented; R bootstrap modes are omitted. |
| `R/predict.R` | `forecast_garch`, implemented through conditional Monte Carlo for all models. |
| `R/persistance.R` | `persistence` and `half_life`. |
| `R/unconditional.R` | `unconditional_variance` and effective-intercept calculations. |
| `R/newsimpact.R` | `news_impact`. Plot construction is omitted. |
| `R/equations.R` | `model_equation`. |
| `R/sandwich.R` | Numerical scores, Hessian, inverse-Hessian, OPG, sandwich covariance, and confidence intervals. HAC/Newey-West covariance is not implemented. |
| `R/profile.R` | `profile_likelihood`. Parallel execution and R summaries are omitted. |
| `R/backtest.R` | `backtest_var` plus Kupiec and Christoffersen coverage tests. Date-index and table output are omitted. |
| `R/benchmark.R` | Exact FCP/Laurent benchmark constants and log-relative-error helper in `tsgarch_benchmarks.f90`. R data loading and formatted tables are omitted. |
| `R/extractors.R` | Replaced by public fields in `garch_fit` and `garch_filter_result`. |
| `R/multispec.R` | Use ordinary arrays of `garch_spec` and user loops. S3 combination and future-based parallelism are omitted. |
| `R/methods.R`, `R/print.R`, `R/reexports.R`, `R/data.R`, `R/tsgarch-package.R` | R classes, methods, formatting, data registration, and reexports omitted. |
| `R/utilities.R` | Relevant specification copying and validation are represented by typed assignment and validation. Date-series conversion helpers are omitted. |
| `RcppExports.R`, `src/RcppExports.cpp`, TMB registration files | R/TMB glue omitted. |
| TMB GARCH/GJR/APARCH/EGARCH/family/component headers | Recursions independently translated into `tsgarch_model.f90`, `tsgarch_simulation.f90`, and related modules. |
| `src/TMB/realgarch.hpp.cpp` | Not included: real-GARCH is not listed by upstream `valid_garch_models()` and has no exported specification path in version 1.0.4. |
