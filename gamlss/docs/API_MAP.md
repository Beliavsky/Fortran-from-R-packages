# API map

| Upstream computational role | Fortran v0.9.0 |
|---|---|
| `gamlss(..., method=RS())` | `fit_gamlss_model(..., GAMLSS_METHOD_RS)` |
| `gamlss(..., method=CG())` | `fit_gamlss_model(..., GAMLSS_METHOD_CG)` |
| `mixed()` | `GAMLSS_METHOD_MIXED` |
| `gamlss.control()` | `gamlss_control_t` |
| parameter fitted values / LP | `result%mu/sigma/nu/tau%fitted`, `%eta` |
| parameter coefficients | `%coefficients` |
| global / penalized deviance | `%global_deviance`, `%penalized_deviance` |
| EDF / df.fit | parameter `%edf`, `result%df_fit` |
| AIC / GAIC / SBC | `result%aic`, `gaic`, `result%sbc` |
| `pb()` / `ps()` core | `fit_p_spline_basis`, `difference_penalty` |
| spline prediction | `predict_p_spline_basis` |
| `cs()` numerical basis role | `natural_spline_design` |
| `cy()` penalty role | `cyclic_difference_penalty` |
| `ridge()` | `ridge_penalty` |
| `fp()` / `gamlss.fp()` | `select_fractional_polynomial`, `predict_fractional_polynomial` |
| `lo()` / `gamlss.lo()` | `fit_loess`, `predict_loess` |
| `pvc()` numerical role | `varying_coefficient_p_spline` |
| `pbm()` numerical role | `fit_monotone_p_spline`, `predict_monotone_p_spline` |
| `random()` / simple `ri()` | `random_intercept_design` |
| `re()` random-intercept role | `fit_gamlss_random_intercept` |
| general `re()` grouped random slopes | `fit_gamlss_random_effects` |
| Gaussian `gls()` residual correlation/variance role | `fit_gamlss_no_gls` |
| random intercepts on multiple GAMLSS parameters | `fit_gamlss_multi_random_intercept` |
| parameter-wise stepwise GAIC/BIC | `stepwise_gaic_parameter` |
| worm-plot numerical summary | `worm_plot_diagnostics` |
| leverage/influence diagnostics | `influence_from_hat` |
| residual Jarque-Bera statistic | `jarque_bera_statistic` |
| `Surv(..., type="interval2")` bridge | `surv_interval2` |
| `Surv(start,stop,event)` bridge | `surv_counting_process` + `entry` |
| delayed-entry censored likelihood | `fit_gamlss_censored(..., entry=...)` |
| additive `pb()` term composition | `build_additive_p_splines` |
| 2D tensor smooth numerical role | `tensor_p_spline_2d` |
| backward/both stepwise GAIC | `stepwise_gaic_mu` |
| case bootstrap | `bootstrap_gamlss_cases` |
| profile LR interval | `profile_likelihood_ci` |
| compiled `genD()` | `all_pair_difference_matrix` |
| `pcat()` / `gamlss.pcat()` | `fit_pcat`, `pcat_fused_groups` |
| censored GAMLSS family fit | `fit_gamlss_censored` |
| `Surv(time,event)` bridge | `surv_right_censoring` |
| hat values | `hat_values_penalized` |
| `rqres()` | `randomized_quantile_residuals` |
| residual ACF | `residual_acf` |
| `devianceIncr()` | `deviance_increment` |
| LR statistic | `likelihood_ratio_stat` |
| `fitDist` / `chooseDist` numerical comparison | `compare_families`, `best_family` |
| stepwise GAIC numerical role | `forward_gaic_mu` |
| profile likelihood | `profile_gamlss_coefficient` |
| `lms()` | `fit_lms` |
| LMS prediction | `predict_lms` |
| centile computation | `lms_centiles` |
| `predictAll()` parameter prediction | `predict_gamlss_parameters` |

Distribution constants and `d/p/q/r` APIs are exported through the umbrella
`use gamlss` from the vendored `gamlss_dist` module.


## v0.5 additions

| Fortran API | Computational role |
|---|---|
| `fit_gamlss_correlated_rs` | Parameter-wise correlated RS/Fisher working-response fitting for non-Gaussian families |
| `correlated_rs_result_t` | GAMLSS result plus fitted/shared `nlme` correlation and base variance scale |
| `fit_gamlss_multi_random_effects` | Simultaneous q-dimensional grouped random-effect blocks on multiple distribution parameters |
| `multi_random_effects_result_t` | Group effects plus per-parameter covariance/precision matrices |
| `cross_validate_gamlss` | Explicit-fold out-of-sample log-likelihood validation |
| `randomized_quantile_residuals_all` | Quantile residuals through the full generic family CDF dispatcher |


## v0.6 additions

| Fortran API | Computational role |
|---|---|
| `fit_gamlss_gaussian_copula` | Exact Gaussian-copula joint likelihood for continuous GAMLSS margins |
| `gaussian_copula_result_t` | Joint/marginal/copula likelihoods, Gaussian scores, fitted correlation and covariance |
| `fit_gamlss_joint_random_effects` | Joint grouped random effects with cross-parameter covariance |
| `joint_random_effects_result_t` | Group effects and full covariance/precision over active parameter x random-term blocks |


## v0.7 additions

| Fortran API | Computational role |
|---|---|
| `mvn_rectangle_probability` | Multivariate-normal rectangle probability by Genz sequential QMC integration |
| `mvn_conditional` | Conditional Gaussian mean/covariance for mixed copula likelihoods |
| `fit_gamlss_gaussian_copula_mixed` | Gaussian-copula joint likelihood for discrete and mixed atomic/continuous margins |
| `gaussian_copula_mixed_result_t` | Joint likelihood, group contributions, fitted correlation and parameter covariance |
| `fit_gamlss_joint_random_effects_ghq` | Direct Gauss-Hermite marginal likelihood for low-dimensional cross-parameter random effects |
| `joint_random_ghq_result_t` | Marginal likelihood, full random covariance and quadrature posterior effects |
| `family_cdf_left` | Left-limit CDF used to construct latent rectangles |
| `family_observation_is_atom` | Detect discrete/endpoint atomic observations |


## v0.8 additions

| Fortran API | Computational role |
|---|---|
| `fit_gamlss_joint_random_effects_ais` | Adaptive QMC importance-sampling marginal likelihood for moderate-dimensional joint random effects |
| `joint_random_ais_result_t` | Marginal likelihood, full random covariance, posterior means/covariances, group ESS diagnostics |


## v0.9 additions

| Fortran API | Computational role |
|---|---|
| `marginal_predict_eta` | Matrix-first `getMarginal()` calculation from a random-effect-free linear predictor |
| `get_marginal_random_intercept` | `getMarginal()` adapter for `random_intercept_result_t` and group labels |
| `MARGINAL_INTEGRATE` | Upstream `method="integrate"` Gaussian random-effect integration |
| `MARGINAL_QFUNCTION` | Upstream 999-normal-quantile deterministic averaging |
| `MARGINAL_RANDOM` | Upstream Monte Carlo averaging, default 10,000 draws |
| `MARGINAL_NONE` | Remove the random term and apply the inverse link only |
