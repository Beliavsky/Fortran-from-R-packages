# API map

This table maps the non-graphical computational surface of `tsgarch` 1.0.4 to
the Fortran library.

| R package surface | Fortran equivalent | Notes |
|---|---|---|
| `garch_modelspec()` | `type(garch_spec)` | Set model, distribution, orders, initialization, targeting, and related options directly. |
| model parameter vectors | `type(garch_parameters)` | Arrays hold ARCH, GARCH, asymmetry, family-GARCH, and regressor coefficients. |
| `estimate.tsgarch.spec()` | `estimate_garch()` | Bounded likelihood optimization and optional numerical inference. |
| `tsfilter.tsgarch.spec()` | `filter_garch()` | Returns conditional sigma, variance, residuals, log-likelihood contributions, and components. |
| `tsfilter.tsgarch.estimate()` | `filter_garch()` using `fit%spec` and `fit%parameters` | No S3 dispatch is required. |
| `simulate.tsgarch.spec()` | `simulate_garch()` | Unconditional simulation with burn-in and multiple paths. |
| conditional simulation used by prediction | `simulate_conditional()` | Starts from observed history. |
| `predict.tsgarch.estimate()` | `forecast_garch()` | Monte Carlo mean, variance, sigma, and requested quantiles. |
| `persistence()` | `persistence()` | Model- and distribution-aware persistence. |
| `omega()` | `effective_omega()` or `effective_sample_omega()` | Accounts for variance targeting and regressors. |
| `unconditional()` | `unconditional_variance()` | Returns the model's unconditional variance when defined. |
| `halflife()` | `half_life()` | Persistence half-life. |
| `newsimpact()` | `news_impact()` | Evaluates the conditional-variance response over supplied shocks. |
| `tsequation()` | `model_equation()` | Returns a concise equation description. |
| `pit()` | `probability_integral_transform()` | Conditional distribution values for fitted residuals. |
| `logLik()` | `fit%log_likelihood` | Scalar field. |
| `AIC()` | `fit%aic` | Scalar field. |
| `BIC()` | `fit%bic` | Scalar field. |
| `coef()` | `fit%packed_parameters`, `fit%parameter_names` | Natural-parameter ordering is documented by the returned names. |
| `sigma()` | `fit%filtered%sigma` | Vector field. |
| `residuals()` | `fit%filtered%residuals` | Standardized residuals are also returned. |
| `fitted()` | `data - fit%filtered%residuals` | The library stores the components needed to construct fitted values. |
| `nobs()` | `fit%filtered%nobs` | Integer field. |
| `estfun()` | `fit%scores` | Per-observation numerical score matrix. |
| `bread()` | `fit%hessian` | Numerical observed Hessian of the negative log likelihood. |
| `vcov()` | `fit%covariance` | Inverse-Hessian covariance when available. |
| OPG covariance | `covariance_opg()` | Uses the observation score matrix. |
| robust sandwich covariance | `covariance_sandwich()` | Hessian bread and OPG meat. |
| `confint()` | `confidence_intervals()` | Wald intervals on packed natural parameters. |
| `tsprofile()` | `profile_likelihood()` | Evaluates likelihood on a user-provided parameter grid. |
| `tsbacktest()` | `backtest_var()` | Rolling or expanding one-step VaR forecasts and coverage tests. |
| Kupiec and Christoffersen tests | `coverage_tests()` | Also called by `backtest_var()`. |
| `benchmark_fcp()` | `fcp_benchmark_data()` | Published coefficient and standard-error fixtures. |
| `benchmark_laurent()` | `laurent_benchmark_data()` | Published coefficient and standard-error fixtures. |
| benchmark relative-error helper | `log_relative_error()` | Elemental scalar function. |
| `nloptr_fast_options()` / `nloptr_global_options()` | `type(fit_options)` | The exact nloptr option lists do not map to the internal optimizer. |
| `to_multi_estimate()` | omitted | R class/list conversion without a numerical kernel. Use arrays of `garch_fit`. |
| multi-spec parallel estimation | user loop over `estimate_garch()` | Parallel scheduling is outside the numerical library. |
| plotting and `as_flextable()` methods | omitted | Presentation-only code. |
| `xts`/`zoo` date handling | omitted | Pass plain ordered numeric arrays; retain dates in the calling application. |

## Model identifiers

Use the following lowercase `spec%model` values:

`garch`, `gjrgarch`, `aparch`, `egarch`, `fgarch`, `cgarch`, `igarch`, `ewma`.

Use these lowercase `spec%distribution` values:

`norm`, `std`, `snorm`, `sstd`, `ged`, `sged`, `nig`, `gh`, `jsu`, `ghst`.
