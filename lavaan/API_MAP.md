# API map

| Upstream lavaan area | Fortran API |
|---|---|
| RAM implied covariance / means | `ram_sigma`, `ram_mu` |
| LISREL/RAM representation | `ram_model`, `ram_from_lisrel` |
| free parameter mapping | `ram_free_map`, `ram_get_free`, `ram_set_free` |
| ML / GLS / ULS / WLS / DWLS objectives | `objective_ml`, `objective_gls`, `objective_uls`, `objective_wls`, `objective_dwls` |
| raw normal likelihood | `mvn_loglik_complete` |
| missing-data likelihood | `mvn_loglik_missing`, `fit_ram_fiml` |
| covariance/raw SEM fitting | `fit_ram_cov`, `fit_ram_data` |
| multigroup equality/free systems | `ram_group_spec`, `independent_group_links`, `fit_ram_multigroup_cov` |
| multigroup raw-data FIML | `ram_group_data`, `fit_ram_multigroup_fiml` |
| nonlinear/equality/inequality constraints | `constraint_callback`, `fit_ram_cov_constrained` |
| EFA extraction | `efa_principal_axis`, `efa_ml`, `efa_fit_cov` |
| EFA rotation | `efa_fit_cov(..., rotation=...)` through vendored GPArotation |
| SAM | `sam_fit_cov`, `sam_fit_data`, `sam_fix_measurement`, `sam_propagate_uncertainty`, `sam_yuan_chan_test`, `sam_continuous_gamma`, `sam_browne_unbiased_gamma` |
| MIIV / IV 2SLS | `miiv_2sls`, `miiv_2sls_cov`, `miiv_estimate_uls`, `miiv_estimate_gls`, `miiv_estimate_2rls`, `miiv_estimate_rls` |
| MIIV covariance linearization | `miiv_jacobian_uls`, `miiv_jacobian_gls`, `miiv_jacobian_2rls`, `miiv_jacobian_rls`, `miiv_vcov_from_gamma` |
| model-implied IV screening/equations | `ram_miiv_candidates`, `ram_miiv_equations`, `ram_miiv_marker_equations`, `miiv_proxy_node` |
| mixed/conditional-x Muthen statistics | `muthen1984_mixed` |
| automatic/name-aware MIIV marker rewriting | `miiv_auto_markers`, `ram_miiv_named_equations`, `miiv_partable_markers` |
| SAM joint block covariance / second-order bias | `sam_block_covariance`, `sam_second_order_bias` |
| missing-response random coefficient ML | `random_coefficient_loglik_missing`, `fit_random_coefficient_missing_ml`, `random_effects_eb_missing` |
| robust nested-model differences | `satorra_bentler_difference_2001`, `satorra_bentler_difference_2010` |
| robust/sandwich SE | `robust_ml_inference` |
| clustered sandwich SE | `robust_ml_inference(..., cluster=...)` |
| Satorra-Bentler covariance correction | `robust_sem_result%sb_scaling`, `%chisq_scaled` |
| scaled robust tests | `covariance_scaled_tests`, `scaled_tests_from_ugamma`, `yuan_bentler_from_traces`, `hayakawa_trace_corrected`, `hayakawa_adjusted_tests`, `browne_residual_test`, `browne_residual_nt` |
| modification indices / score test | `modification_indices_cov` |
| Wald tests | `wald_test` |
| model standardization | `standardized_ram` |
| covariance residuals | `residual_covariance` |
| regression factor scores | `factor_scores_regression` |
| ordinal thresholds | `ordinal_thresholds` |
| bivariate ordinal probabilities | `bvn_rectangle` |
| polychoric correlation | `polychoric_table`, `polychoric_matrix` |
| correlation-only categorical WLS | `ordinal_wls_correlation_weights` |
| threshold + correlation categorical WLS | `categorical_wls_statistics`, `categorical_wls_statistics_analytic`, `muthen1984_ordinal`, `categorical_wls_statistics_muthen` |
| ordinal pairwise likelihood | `fit_ram_pml_ordinal` |
| mixed continuous/ordinal PML | `fit_ram_pml_mixed` |
| ordinal factor MML | `fit_mml_ordinal_factor`, `mml_ordinal_loglik` |
| mixed continuous/ordinal MML | `fit_mml_mixed_factor`, `mml_mixed_loglik` |
| high-dimensional deterministic MML | `fit_mml_mixed_factor_qmc`, `mml_mixed_loglik_qmc`, `halton_normal_nodes` |
| posterior-adaptive MML | `fit_mml_mixed_factor_adaptive`, `mml_mixed_loglik_adaptive` |
| Gauss-Hermite nodes/weights | `gauss_hermite_normal` |
| bootstrap inference | `bootstrap_ram_data` |
| clustered/two-level moment SEM | `fit_ram_twolevel` |
| exact complete-data two-level ML | `fit_ram_twolevel_ml` |
| missing-data two-level FIML + H1 | `fit_ram_twolevel_fiml` |
| random-intercept/random-slope Gaussian ML | `fit_random_coefficient_ml`, `random_coefficient_loglik`, `random_effects_eb` |
| RAM simulation | `simulate_ram` |
| `lav_matrix_vec` / `vech` | `vec`, `vech`, `vech_reverse` |
| duplication/commutation matrices | `duplication_matrix`, `commutation_matrix` |
| symmetric square roots | `symmetric_sqrt` |
| orthogonal complement | `orthogonal_complement` |

R syntax parsing, S4 dispatch, plotting, printing, and data-frame/model-frame plumbing are intentionally outside the Fortran API.
