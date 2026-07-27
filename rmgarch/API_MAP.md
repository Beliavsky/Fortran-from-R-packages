# API map

This map describes only procedures that are present in this release. "Tested"
means an included test program calls the path and checks numerical invariants or
finite output. "Partial" means the numerical idea is available but the full R
options and object behavior are not.

| R API/concept | Fortran API | Verified status |
|---|---|---|
| `dccspec` | `type(dcc_spec)`, `make_dcc_spec` | Tested; numerical specification only |
| `dccfit` | `fit_dcc`, `fit_dcc11`, `fit_two_step_dcc_general` | Tested; scalar DCC/ADCC with Gaussian GARCH(1,1) margins |
| `dccfilter` | `dcc_filter` | Tested |
| `dccforecast` | `dcc_forecast`, `dcc_forecast_history` | Tested |
| `dccsim` | `simulate_dcc` | Tested for Gaussian, Student, and Laplace |
| `dccroll` | `roll_dcc`, `roll_two_step_dcc` | Tested sequential rolling core; no parallel/refit persistence |
| FDCC fit/filter/forecast/sim | `fit_fdcc11`, `fdcc_filter`, `fdcc_forecast`, `simulate_fdcc` | Tested grouped `(1,1)` implementation |
| `cgarchfit` | `fit_copula`, `fit_copula_garch11` | Tested Gaussian model-level workflow; Gaussian/Student copula primitives tested |
| `cgarchfilter` | `filter_copula_garch11` | Tested for dynamic Gaussian copula-GARCH |
| `cgarchsim` | `simulate_fitted_copula_garch11` | Tested for dynamic Gaussian copula-GARCH |
| Copula transformations | `parametric_uniform_transform`, `empirical_uniform_transform`, `copula_score_transform`, `score_to_uniform_transform` | Tested |
| `gogarchfit` | `fit_gogarch11` | Tested square ICA plus Gaussian GARCH(1,1) factors |
| `gogarchfilter` | `filter_gogarch11` | Tested |
| `gogarchforecast` | `forecast_gogarch11` | Tested |
| `gogarchsim` | `simulate_gogarch`, `simulate_fitted_gogarch11` | Tested |
| `gogarchroll` | `roll_gogarch11` | Tested sequential rolling core |
| `rcov`, `rcor`, `sigma` | DCC/GO-GARCH result arrays; `gogarch_covariance`, `gogarch_correlation`, `gogarch_sigma` | Tested |
| `rcoskew`, `rcokurt` | `gogarch_coskewness`, `gogarch_cokurtosis`, `gogarch_moments_at` | Tested |
| `gportmoments`, `nportmoments` | `portfolio_factor_moments` | Tested moment contraction core |
| `betacovar`, `betacoskew`, `betacokurt` | `portfolio_covariance_beta`, `portfolio_coskew_beta`, `portfolio_cokurt_beta` | Tested |
| `fastica` | `fastica` | Tested reconstruction |
| `radical` | `radical` | Tested reconstruction |
| `wmargin` | `weighted_margin`, `weighted_margin_path` | Tested for Student and Laplace paths; Normal shares the same contraction core |
| `dfft`, `pfft`, `qfft`, convolution | `fft_transform`, `convolve_grid_distributions`, `grid_density`, `grid_cdf`, `grid_quantile`, `grid_moments` | Tested through Normal-grid convolution |
| `fmoments`, `fscenario` | `scenario_moments`, `simulate_dcc_scenarios`, `portfolio_scenarios` | Tested DCC scenario/moment core; not full R class behavior |
| `varxfit` | `fit_varx`, `fit_varx_robust` | Tested |
| `varxfilter` | `filter_varx` | Tested |
| `varxforecast` | `forecast_varx` | Tested |
| `varxsim` | `simulate_varx` | Tested |
| `DCCtest` | `dcc_constancy_test` | Tested |
| Mardia diagnostics | `mardia_test` | Tested |
| `cordist` | `correlation_distance`, `correlation_distance_matrix` | Present; distance core covered by compilation, plotting omitted |
| Multivariate distributions | `multivariate_*_logpdf`, `random_multivariate_*` | Tested for Normal, Student, Laplace |
| `rshape`, `rskew` | Stored scalar shape in fit/spec types; no R-style accessor wrapper | No class wrapper needed; skewed DCC family not implemented |
| Plotting and R methods | - | Intentionally omitted |
