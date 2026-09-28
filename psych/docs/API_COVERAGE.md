# API coverage

Coverage basis: exported computational R functions, following the translation conventions in the repository.

- Package status: **substantial**
- Exported computational functions: **344**
- Mapped functions: **344**
- Mapping fraction: **100.0%**
- Unmapped computational functions: **0**

The denominator excludes plotting, printing, formatting, dataset/data-download, interactive-display, presentation-only exports, R-only syntax/object constructors, pure selection/reshaping helpers, and functions that merely dispatch computation without adding a numerical operation. A mapped function may still be partial or substantial rather than R-interface complete.

A coverage re-audit removed ten previously counted non-computational exports: `acs`, `cs`, `errorCircles`, `fromTo`, `histBy`, `irt.select`, `make.irt.stats`, `mlArrange`, `phi.list`, and `structure.sem`.

## Untranslated computational exports

None in the current 344-function computational coverage basis. A mapped function can still be partial or substantial rather than R-interface complete.

## Mapped functions

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `%+%` | `psych_matrix` | `matrix_addition_operator` | substantial |
| `alpha` | `psych_reliability` | `cronbach_alpha_data` | partial |
| `alpha.ci` | `psych_inference_extended` | `alpha_confidence_interval` | substantial |
| `alpha2r` | `psych_basic` | `alpha_to_r` | complete |
| `AUC` | `psych_basic` | `auc_from_table` | substantial |
| `autoR` | `psych_correlations` | `auto_correlation` | substantial |
| `bassAckward` | `psych_more_extended` | `bass_ackward_factor` | substantial |
| `bestItems` | `psych_selection_extended` | `best_items_from_correlation` | substantial |
| `bestScales` | `psych_more_extended` | `best_scales_from_correlation` | partial |
| `bifactor` | `psych_factor_extended` | `bifactor_rotate` | substantial |
| `bigCor` | `psych_data_extended` | `big_correlation_matrix` | substantial |
| `biquartimin` | `psych_factor_extended` | `biquartimin_rotate` | substantial |
| `biserial` | `psych_correlations` | `biserial_correlation` | substantial |
| `block.random` | `psych_additional` | `block_random_design` | substantial |
| `cd.validity` | `psych_validity_extended` | `effect_validity_by_keys` | substantial |
| `chi2r` | `psych_basic` | `chi2_to_r` | complete |
| `circ.sim` | `psych_irt_extended` | `simulate_item_circumplex` | substantial |
| `circ.simulation` | `psych_selection_extended` | `simulate_circumplex_benchmark` | substantial |
| `circ.tests` | `psych_circular` | `circ_tests_values` | substantial |
| `circadian.cor` | `psych_circular` | `circadian_correlation_matrix` | substantial |
| `circadian.F` | `psych_circular` | `circadian_f_test` | partial |
| `circadian.linear.cor` | `psych_circular` | `circadian_linear_correlation` | substantial |
| `circadian.mean` | `psych_circular` | `circular_mean_value` | substantial |
| `circadian.phase` | `psych_circular` | `cosinor_fit` | partial |
| `circadian.reliability` | `psych_computational_extended` | `circadian_reliability_split` | substantial |
| `circadian.sd` | `psych_circular` | `circadian_stats_vector` | substantial |
| `circadian.stats` | `psych_circular` | `circadian_stats_vector` | substantial |
| `circular.cor` | `psych_circular` | `circular_correlation_matrix` | substantial |
| `circular.mean` | `psych_circular` | `circular_mean_value` | complete |
| `cluster.cor` | `psych_scoring_extended` | `cluster_correlation_summary` | substantial |
| `cluster.fit` | `psych_additional` | `cluster_fit_metrics` | substantial |
| `cluster.loadings` | `psych_scoring_extended` | `cluster_loadings_summary` | substantial |
| `cluster2keys` | `psych_scoring_extended` | `cluster_assignments_to_keys` | complete |
| `cohen.d` | `psych_correlations` | `cohen_d_value` | partial |
| `cohen.d.by` | `psych_selection_extended` | `cohen_d_by_groups` | partial |
| `cohen.d.ci` | `psych_inference_extended` | `effect_size_confidence_interval` | partial |
| `cohen.kappa` | `psych_correlations` | `cohen_kappa_value` | substantial |
| `cohen.profile` | `psych_additional` | `cohen_profile_similarity` | substantial |
| `cohenPooled` | `psych_correlations` | `pooled_sd_value` | substantial |
| `comorbidity` | `psych_additional` | `comorbidity_statistics` | substantial |
| `congeneric.sim` | `psych_simulation` | `congeneric_correlation`, `simulate_correlated_normal` | partial |
| `congruence` | `psych_correlations` | `congruence_matrix` | substantial |
| `cor.ci` | `psych_basic` | `correlation_ci` | partial |
| `cor.smooth` | `psych_correlations` | `smooth_correlation` | substantial |
| `cor.smoother` | `psych_correlations` | `smooth_correlation` | substantial |
| `cor.wt` | `psych_correlations` | `weighted_correlation_matrix` | substantial |
| `cor2` | `psych_additional` | `cor2_matrix` | partial |
| `cor2cov` | `psych_basic` | `r_to_covariance` | complete |
| `cor2dist` | `psych_correlations` | `correlation_to_distance` | complete |
| `corCi` | `psych_basic` | `correlation_ci` | partial |
| `corFiml` | `psych_more_extended` | `correlation_fiml` | substantial |
| `corInfo` | `psych_inference_extended` | `correlation_information` | substantial |
| `corr.p` | `psych_correlations` | `correlation_significance` | substantial |
| `corr.test` | `psych_correlations` | `correlation_significance` | partial |
| `correct.cor` | `psych_correlations` | `correct_correlation` | substantial |
| `corTest` | `psych_correlations` | `correlation_significance` | partial |
| `cortest` | `psych_inference_extended` | `correlation_normal_test` | partial |
| `cortest.bartlett` | `psych_correlations` | `bartlett_correlation_test` | substantial |
| `cortest.jennrich` | `psych_inference_extended` | `jennrich_correlation_test` | substantial |
| `cortest.mat` | `psych_inference_extended` | `correlation_matrix_test` | substantial |
| `cortest.normal` | `psych_inference_extended` | `correlation_normal_test` | substantial |
| `cosinor` | `psych_circular` | `cosinor_fit` | partial |
| `cosinor.period` | `psych_circular` | `cosinor_period_fit` | partial |
| `count.pairwise` | `psych_descriptive` | `pairwise_count_matrix` | substantial |
| `crossValidation` | `psych_validity_extended` | `cross_validation_from_correlation` | partial |
| `crossValidationBoot` | `psych_computational_extended` | `cross_validation_bootstrap` | partial |
| `cta` | `psych_more_extended` | `cta_dynamics` | substantial |
| `cta.15` | `psych_more_extended` | `cta15_dynamics` | substantial |
| `d.ci` | `psych_inference_extended` | `effect_size_confidence_interval` | partial |
| `d.robust` | `psych_inference_extended` | `robust_effect_size_two_groups` | substantial |
| `d2CL` | `psych_basic` | `d_to_common_language` | complete |
| `d2OVL` | `psych_basic` | `d_to_overlap` | complete |
| `d2OVL2` | `psych_basic` | `d_to_overlap2` | complete |
| `d2r` | `psych_basic` | `d_to_r` | complete |
| `d2t` | `psych_basic` | `d_to_t` | substantial |
| `d2U3` | `psych_basic` | `d_to_u3` | complete |
| `describe` | `psych_descriptive` | `describe_vector`, `describe_matrix` | substantial |
| `describe.by` | `psych_scoring_extended` | `grouped_describe_matrix` | partial |
| `describeBy` | `psych_scoring_extended` | `grouped_describe_matrix` | partial |
| `describeData` | `psych_data_extended` | `describe_data_numeric` | partial |
| `describeFast` | `psych_data_extended` | `describe_fast_numeric` | partial |
| `directSl` | `psych_omega_extended` | `direct_schmid_leiman` | substantial |
| `distance` | `psych_correlations` | `profile_distance_matrix` | substantial |
| `dummy.code` | `psych_scoring` | `dummy_code_matrix` | substantial |
| `eigen.loadings` | `psych_factor` | `eigen_loadings` | complete |
| `eigenCi` | `psych_additional` | `eigen_confidence_intervals` | substantial |
| `equamax` | `psych_factor_extended` | `equamax_rotate` | substantial |
| `esem` | `psych_more_extended` | `exploratory_sem` | substantial |
| `fa` | `psych_factor` | `principal_axis_factor` | partial |
| `fa.congruence` | `psych_correlations` | `congruence_matrix` | substantial |
| `fa.extend` | `psych_factor_extended` | `factor_extension_loadings` | substantial |
| `fa.extension` | `psych_factor_extended` | `factor_extension_loadings` | substantial |
| `fa.multi` | `psych_computational_extended` | `factor_analysis_multi` | substantial |
| `fa.organize` | `psych_computational_extended` | `organize_factor_loadings` | substantial |
| `fa.parallel` | `psych_parallel_extended` | `parallel_analysis_normal` | partial |
| `fa.parallel.poly` | `psych_more_extended` | `parallel_analysis_poly` | substantial |
| `fa.poly` | `psych_computational_extended` | `factor_analysis_poly_data` | partial |
| `fa.pooled` | `psych_more_extended` | `factor_analysis_pooled` | substantial |
| `fa.random` | `psych_parallel_extended` | `random_effect_factor` | substantial |
| `fa.sapa` | `psych_more_extended` | `factor_analysis_sapa` | partial |
| `fa.sort` | `psych_factor_extended` | `sort_factor_loadings` | substantial |
| `fa.stats` | `psych_factor_extended` | `factor_fit_statistics` | partial |
| `fa2irt` | `psych_factor_extended` | `factor_to_irt_parameters` | substantial |
| `faBy` | `psych_more_extended` | `factor_analysis_by_groups` | substantial |
| `fac` | `psych_factor_extended` | `minimum_residual_factor` | partial |
| `faCor` | `psych_computational_extended` | `factor_solution_correlation` | partial |
| `factor.congruence` | `psych_correlations` | `congruence_matrix` | substantial |
| `factor.fit` | `psych_correlations` | `factor_fit_value` | partial |
| `factor.minres` | `psych_factor_extended` | `minimum_residual_factor` | substantial |
| `factor.model` | `psych_correlations` | `factor_model_matrix` | substantial |
| `factor.pa` | `psych_factor` | `principal_axis_factor` | substantial |
| `factor.residuals` | `psych_factor` | `factor_residuals` | substantial |
| `factor.rotate` | `psych_factor` | `varimax_rotate`, `promax_rotate`, `target_rotation` | partial |
| `factor.scores` | `psych_factor` | `factor_scores_regression` | partial |
| `factor.stats` | `psych_factor_extended` | `factor_fit_statistics` | partial |
| `factor.wls` | `psych_factor_extended` | `weighted_least_squares_factor` | substantial |
| `factor2cluster` | `psych_factor` | `factor_to_clusters` | substantial |
| `faReg` | `psych_factor_extended` | `factor_regression_partition` | partial |
| `faRegression` | `psych_factor_extended` | `factor_regression_partition` | partial |
| `faRotate` | `psych_factor` | `varimax_rotate`, `promax_rotate`, `target_rotation` | partial |
| `faRotations` | `psych_factor` | `varimax_rotate`, `promax_rotate`, `target_rotation` | partial |
| `fisherz` | `psych_basic` | `fisher_z` | complete |
| `fisherz2r` | `psych_basic` | `fisher_r` | complete |
| `fsi` | `psych_additional` | `factor_score_indeterminacy` | substantial |
| `g2r` | `psych_basic` | `g_to_r` | complete |
| `geometric.mean` | `psych_basic` | `geometric_mean` | substantial |
| `glb` | `psych_omega_extended` | `glb_reliability_summary` | partial |
| `glb.algebraic` | `psych_more_extended` | `glb_algebraic_admm` | substantial |
| `glb.fa` | `psych_omega_extended` | `glb_factor_reliability` | substantial |
| `guttman` | `psych_reliability` | `guttman_lambdas` | substantial |
| `harmonic.mean` | `psych_basic` | `harmonic_mean` | substantial |
| `ICC` | `psych_reliability_extended` | `icc_complete_matrix` | partial |
| `ICLUST` | `psych_iclust_extended` | `iclust_hierarchical` | partial |
| `iclust` | `psych_iclust_extended` | `iclust_hierarchical` | partial |
| `ICLUST.cluster` | `psych_iclust_extended` | `iclust_hierarchical` | substantial |
| `ICLUST.sort` | `psych_iclust_extended` | `iclust_sort_loadings` | substantial |
| `iclust.sort` | `psych_iclust_extended` | `iclust_sort_loadings` | substantial |
| `interbattery` | `psych_additional` | `interbattery_loadings` | partial |
| `interp.median` | `psych_interpolation` | `interpolated_quantile` | substantial |
| `interp.q` | `psych_interpolation` | `interpolated_quantile` | partial |
| `interp.quantiles` | `psych_interpolation` | `interpolated_quantile` | partial |
| `interp.quart` | `psych_interpolation` | `interpolated_quartiles` | partial |
| `interp.quartiles` | `psych_interpolation` | `interpolated_quartiles` | partial |
| `interp.values` | `psych_interpolation` | `interpolated_values` | substantial |
| `irt.0p` | `psych_irt` | `irt_zero_parameter` | substantial |
| `irt.1p` | `psych_irt` | `rasch_person_scores` | substantial |
| `irt.2p` | `psych_irt` | `two_parameter_person_scores` | substantial |
| `irt.discrim` | `psych_irt` | `irt_discrimination` | substantial |
| `irt.fa` | `psych_computational_extended` | `irt_factor_analysis` | substantial |
| `irt.item.diff.rasch` | `psych_irt` | `rasch_item_difficulty` | substantial |
| `irt.person.rasch` | `psych_irt` | `rasch_person_scores` | substantial |
| `irt.responses` | `psych_irt_extended` | `irt_response_curves` | partial |
| `irt.se` | `psych_irt_extended` | `irt_standard_errors` | substantial |
| `irt.stats.like` | `psych_irt_extended` | `irt_stats_from_loadings` | partial |
| `irt.tau` | `psych_irt_extended` | `irt_thresholds` | substantial |
| `isCorrelation` | `psych_correlations` | `is_correlation_matrix` | substantial |
| `isCovariance` | `psych_correlations` | `is_covariance_matrix` | substantial |
| `item.dichot` | `psych_simulation` | `simulate_item_dichot` | substantial |
| `item.sim` | `psych_irt_extended` | `simulate_item_circumplex` | partial |
| `item.validity` | `psych_validity_extended` | `item_validity_by_keys` | substantial |
| `itemSort` | `psych_scoring_extended` | `item_sort_indices` | partial |
| `kaiser` | `psych_factor` | `kaiser_normalize` | substantial |
| `keysort` | `psych_scoring_extended` | `key_sort_columns` | substantial |
| `KMO` | `psych_correlations` | `kmo_statistic` | substantial |
| `kurtosi` | `psych_basic` | `excess_kurtosis` | substantial |
| `lmCor` | `psych_validity_extended` | `regression_from_matrices` | partial |
| `logistic` | `psych_basic` | `logistic_value` | substantial |
| `logistic.grm` | `psych_basic` | `graded_response_probability` | substantial |
| `logit` | `psych_basic` | `logit_value` | complete |
| `lowerUpper` | `psych_matrix` | `lower_upper_merge` | substantial |
| `m2d` | `psych_basic` | `means_to_d` | substantial |
| `m2t` | `psych_basic` | `means_to_t` | substantial |
| `make.congeneric` | `psych_simulation` | `make_congeneric_model` | substantial |
| `make.hierarchical` | `psych_simulation` | `make_hierarchical_model` | substantial |
| `mardia` | `psych_correlations` | `mardia_statistics` | substantial |
| `mat.regress` | `psych_validity_extended` | `regression_from_matrices` | partial |
| `mat.sort` | `psych_matrix` | `sort_matrix_by_loadings` | substantial |
| `matMult` | `psych_matrix` | `missing_matmul` | substantial |
| `matReg` | `psych_validity_extended` | `regression_from_matrices` | partial |
| `matSort` | `psych_matrix` | `sort_matrix_by_loadings` | substantial |
| `mediate` | `psych_more_extended` | `mediate_from_covariance` | partial |
| `minkowski` | `psych_basic` | `minkowski_distance` | substantial |
| `mixed.cor` | `psych_data_extended` | `mixed_correlation_matrix` | substantial |
| `mixedCor` | `psych_data_extended` | `mixed_correlation_matrix` | substantial |
| `mlr` | `psych_computational_extended` | `multilevel_reliability_balanced` | partial |
| `mssd` | `psych_basic` | `mssd_value` | substantial |
| `multilevel.reliability` | `psych_computational_extended` | `multilevel_reliability_balanced` | partial |
| `nfactors` | `psych_validity_extended` | `vss_complexity_fit` | partial |
| `omega` | `psych_omega_extended` | `omega_from_correlation` | partial |
| `omegaDirect` | `psych_omega_extended` | `omega_direct_from_correlation` | substantial |
| `omegaFromSem` | `psych_computational_extended` | `omega_from_cfa_loadings` | substantial |
| `omegaSem` | `psych_omega_sem` | `omega_sem_bifactor` | substantial |
| `omegah` | `psych_omega_extended` | `omega_from_correlation` | substantial |
| `omegaStats` | `psych_reliability_extended` | `omega_statistics` | partial |
| `outlier` | `psych_correlations` | `mahalanobis_squared` | partial |
| `p.rep` | `psych_additional` | `replication_probability` | substantial |
| `p.rep.f` | `psych_additional` | `replication_probability_f` | substantial |
| `p.rep.r` | `psych_additional` | `replication_probability_r` | substantial |
| `p.rep.t` | `psych_additional` | `replication_probability_t` | substantial |
| `paired.r` | `psych_correlations` | `paired_correlation_test` | substantial |
| `pairwiseCount` | `psych_descriptive` | `pairwise_count_matrix` | substantial |
| `pairwiseCountBig` | `psych_scoring_extended` | `pairwise_count_blocked` | substantial |
| `pairwiseDescribe` | `psych_scoring_extended` | `pairwise_describe_counts` | substantial |
| `pairwiseImpute` | `psych_scoring_extended` | `pairwise_impute_by_keys` | substantial |
| `pairwiseSample` | `psych_additional` | `pairwise_sample_description` | substantial |
| `pairwiseZero` | `psych_descriptive` | `pairwise_zero_matrix` | partial |
| `parcels` | `psych_selection_extended` | `parcel_keys_from_correlation` | substantial |
| `partial.r` | `psych_correlations` | `partial_correlation` | substantial |
| `paSelect` | `psych_computational_extended` | `parallel_select_by_keys` | partial |
| `pca` | `psych_factor` | `principal_components_data` | partial |
| `phi` | `psych_basic` | `phi_coefficient` | substantial |
| `phi.demo` | `psych_more_extended` | `phi_demo_simulation` | substantial |
| `phi2poly` | `psych_polychoric` | `phi_to_tetrachoric` | partial |
| `phi2poly.matrix` | `psych_validity_extended` | `phi_to_tetrachoric_values` | substantial |
| `phi2tetra` | `psych_polychoric` | `phi_to_tetrachoric` | partial |
| `Pinv` | `psych_correlations` | `pseudoinverse` | substantial |
| `poly.mat` | `psych_polychoric` | `polychoric_categories` | partial |
| `polychoric` | `psych_polychoric` | `polychoric_categories` | partial |
| `polydi` | `psych_data_extended` | `polydi_correlations` | substantial |
| `polyserial` | `psych_remaining` | `polyserial_correlations` | substantial |
| `predicted.validity` | `psych_validity_extended` | `predicted_validity_summary` | partial |
| `principal` | `psych_factor` | `principal_components_data` | partial |
| `Procrustes` | `psych_factor` | `procrustes_rotate` | substantial |
| `Promax` | `psych_factor` | `promax_rotate` | substantial |
| `r.con` | `psych_basic` | `correlation_ci` | substantial |
| `r.test` | `psych_inference_extended` | `correlation_test_general` | substantial |
| `r2alpha` | `psych_basic` | `r_to_alpha` | complete |
| `r2c` | `psych_basic` | `r_to_covariance` | complete |
| `r2chi` | `psych_basic` | `r_to_chi2` | complete |
| `r2d` | `psych_basic` | `r_to_d` | complete |
| `r2p` | `psych_basic` | `correlation_p_value` | substantial |
| `r2t` | `psych_basic` | `r_to_t` | complete |
| `rangeCorrection` | `psych_correlations` | `range_correction` | substantial |
| `reflect` | `psych_factor` | `reflect_loadings` | partial |
| `reliability` | `psych_reliability` | `reliability_summary` | partial |
| `rescale` | `psych_scoring` | `rescale_matrix` | substantial |
| `response.frequencies` | `psych_descriptive` | `response_frequencies` | substantial |
| `responseFrequency` | `psych_descriptive` | `response_frequencies` | substantial |
| `reverse.code` | `psych_scoring` | `reverse_code_matrix` | substantial |
| `RMSEA` | `psych_basic` | `rmsea_value` | substantial |
| `rmssd` | `psych_basic` | `rmssd_value` | substantial |
| `RV` | `psych_correlations` | `rv_coefficient` | substantial |
| `scaling.fits` | `psych_additional` | `scaling_fit` | substantial |
| `schmid` | `psych_omega_extended` | `schmid_leiman_decomposition` | substantial |
| `score.alpha` | `psych_reliability` | `cronbach_alpha_data` | partial |
| `score.irt` | `psych_irt` | `two_parameter_person_scores` | partial |
| `score.irt.2` | `psych_irt` | `two_parameter_person_scores` | partial |
| `score.irt.poly` | `psych_irt_extended` | `irt_score_poly` | substantial |
| `score.items` | `psych_scoring_extended` | `score_items_summary` | partial |
| `score.multiple.choice` | `psych_scoring` | `score_multiple_choice_matrix` | partial |
| `scoreBy` | `psych_more_extended` | `score_by_correlations` | substantial |
| `scoreFast` | `psych_scoring` | `score_key_matrix` | partial |
| `scoreIrt` | `psych_irt` | `two_parameter_person_scores` | partial |
| `scoreIrt.1pl` | `psych_irt` | `rasch_person_scores` | partial |
| `scoreIrt.2pl` | `psych_irt` | `two_parameter_person_scores` | partial |
| `scoreItems` | `psych_scoring_extended` | `score_items_summary` | partial |
| `scoreOverlap` | `psych_scoring_extended` | `score_overlap_summary` | substantial |
| `scoreVeryFast` | `psych_scoring` | `score_key_matrix` | partial |
| `scoreWtd` | `psych_scoring` | `score_weighted_matrix` | partial |
| `SD` | `psych_basic` | `sd_value` | substantial |
| `set.cor` | `psych_validity_extended` | `regression_from_matrices` | partial |
| `setCor` | `psych_validity_extended` | `regression_from_matrices` | partial |
| `shannon` | `psych_basic` | `shannon_entropy` | substantial |
| `sim` | `psych_simulation_algorithms` | `simulate_structural_model` | substantial |
| `sim.anova` | `psych_simulation_algorithms` | `simulate_anova_data` | partial |
| `sim.bonds` | `psych_remaining` | `simulate_bonds_model` | substantial |
| `sim.circ` | `psych_irt_extended` | `simulate_item_circumplex` | substantial |
| `sim.congeneric` | `psych_simulation` | `congeneric_correlation`, `simulate_correlated_normal` | partial |
| `sim.correlation` | `psych_simulation` | `simulate_correlated_normal` | partial |
| `sim.dichot` | `psych_simulation` | `simulate_dichotomous` | partial |
| `sim.general` | `psych_remaining` | `simulate_general_structure` | substantial |
| `sim.hierarchical` | `psych_simulation_algorithms` | `simulate_hierarchical_data` | substantial |
| `sim.irt` | `psych_simulation` | `simulate_npl_items`, `simulate_npn_items` | substantial |
| `sim.item` | `psych_irt_extended` | `simulate_item_circumplex` | partial |
| `sim.minor` | `psych_simulation_algorithms` | `simulate_minor_factors` | substantial |
| `sim.multi` | `psych_more_extended` | `simulate_multi_subject` | substantial |
| `sim.multilevel` | `psych_simulation_algorithms` | `simulate_multilevel_data` | substantial |
| `sim.npl` | `psych_simulation` | `simulate_npl_items` | substantial |
| `sim.npn` | `psych_simulation` | `simulate_npn_items` | substantial |
| `sim.omega` | `psych_more_extended` | `simulate_omega_study` | substantial |
| `sim.parallel` | `psych_selection_extended` | `simulate_parallel_recovery` | substantial |
| `sim.poly` | `psych_simulation` | `simulate_npl_items`, `simulate_npn_items` | substantial |
| `sim.poly.ideal` | `psych_simulation_algorithms` | `simulate_poly_ideal` | substantial |
| `sim.poly.ideal.npl` | `psych_simulation_algorithms` | `simulate_poly_ideal_logistic` | substantial |
| `sim.poly.ideal.npn` | `psych_simulation_algorithms` | `simulate_poly_ideal_normal` | partial |
| `sim.poly.mat` | `psych_simulation` | `simulate_polychoric_marginals` | substantial |
| `sim.poly.npl` | `psych_simulation` | `simulate_npl_items` | substantial |
| `sim.poly.npn` | `psych_simulation` | `simulate_npn_items` | substantial |
| `sim.rasch` | `psych_simulation` | `simulate_rasch` | substantial |
| `sim.simplex` | `psych_simulation` | `simplex_correlation`, `simulate_correlated_normal` | partial |
| `sim.spherical` | `psych_simulation_algorithms` | `simulate_spherical_items` | substantial |
| `sim.structural` | `psych_simulation` | `structure_correlation`, `simulate_correlated_normal` | partial |
| `sim.structure` | `psych_simulation` | `structure_correlation`, `simulate_correlated_normal` | partial |
| `sim.VSS` | `psych_remaining` | `simulate_vss_data` | substantial |
| `simCor` | `psych_simulation` | `simulate_correlated_normal` | partial |
| `simulation.circ` | `psych_selection_extended` | `simulate_circumplex_benchmark` | substantial |
| `skew` | `psych_basic` | `skewness` | substantial |
| `smc` | `psych_correlations` | `squared_multiple_correlations` | substantial |
| `splitHalf` | `psych_reliability` | `split_half_extremes` | partial |
| `statsBy` | `psych_data_extended` | `stats_by_numeric` | partial |
| `statsBy.boot` | `psych_more_extended` | `stats_by_bootstrap` | substantial |
| `super.matrix` | `psych_matrix` | `super_matrix_pair` | partial |
| `superCor` | `psych_additional` | `super_correlation_matrix` | substantial |
| `superMatrix` | `psych_matrix` | `super_matrix_pair` | partial |
| `t2d` | `psych_basic` | `t_to_d` | substantial |
| `t2r` | `psych_basic` | `t_to_r` | complete |
| `tableF` | `psych_polychoric` | `integer_contingency_table` | substantial |
| `target.rot` | `psych_factor` | `target_rotation` | substantial |
| `TargetQ` | `psych_factor` | `target_rotation` | partial |
| `TargetT` | `psych_factor_extended` | `target_orthogonal_rotate` | substantial |
| `tenberge` | `psych_remaining` | `ten_berge_bounds` | substantial |
| `test.irt` | `psych_more_extended` | `test_irt_benchmark` | partial |
| `testRetest` | `psych_omega_extended` | `test_retest_summary` | partial |
| `tetrachoric` | `psych_polychoric` | `tetrachoric_table` | partial |
| `thurstone` | `psych_remaining` | `thurstone_scaling` | substantial |
| `tr` | `psych_basic` | `trace_value` | complete |
| `unidim` | `psych_validity_extended` | `unidimensionality_statistics` | partial |
| `validityItem` | `psych_validity_extended` | `item_validity_by_keys` | substantial |
| `varimin` | `psych_factor_extended` | `varimin_rotate` | substantial |
| `vech` | `psych_matrix` | `vech_values` | substantial |
| `vgQ.bimin` | `psych_factor_extended` | `bimin_objective_gradient` | complete |
| `vgQ.targetQ` | `psych_factor_extended` | `target_objective_gradient` | substantial |
| `vgQ.varimin` | `psych_factor_extended` | `varimin_objective_gradient` | complete |
| `VSS` | `psych_validity_extended` | `vss_complexity_fit` | partial |
| `vss` | `psych_validity_extended` | `vss_complexity_fit` | partial |
| `VSS.parallel` | `psych_selection_extended` | `vss_parallel_random` | partial |
| `VSS.sim` | `psych_remaining` | `simulate_vss_data` | substantial |
| `VSS.simulate` | `psych_remaining` | `simulate_vss_data` | substantial |
| `vssSelect` | `psych_computational_extended` | `vss_select_by_keys` | partial |
| `winsor` | `psych_basic` | `winsorize`, `winsor_mean` | partial |
| `winsor.mean` | `psych_basic` | `winsorize`, `winsor_mean` | partial |
| `winsor.means` | `psych_basic` | `winsorize`, `winsor_mean` | partial |
| `winsor.sd` | `psych_basic` | `winsor_sd` | substantial |
| `winsor.var` | `psych_basic` | `winsor_variance` | substantial |
| `wkappa` | `psych_correlations` | `weighted_kappa_value` | substantial |
| `Yule` | `psych_basic` | `yule_q`, `yule_y` | substantial |
| `Yule.inv` | `psych_basic` | `yule_inverse` | substantial |
| `Yule2phi` | `psych_basic` | `yule_to_phi` | substantial |
| `Yule2phi.matrix` | `psych_validity_extended` | `yule_to_phi_values` | substantial |
| `Yule2poly` | `psych_validity_extended` | `yule_to_tetrachoric_value` | substantial |
| `Yule2poly.matrix` | `psych_validity_extended` | `yule_to_tetrachoric_values` | substantial |
| `Yule2tetra` | `psych_validity_extended` | `yule_to_tetrachoric_value` | substantial |
| `YuleBonett` | `psych_validity_extended` | `yule_bonett_interval` | substantial |
| `YuleCor` | `psych_validity_extended` | `yule_correlation_matrix` | partial |
