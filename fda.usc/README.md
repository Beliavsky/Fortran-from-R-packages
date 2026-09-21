# fda.usc - modern Fortran translation

This directory translates computational routines from the R package **fda.usc 2.2.0** into modern free-form Fortran for use with FPM. The upstream package is by Manuel Febrero Bande and Manuel Oviedo de la Fuente, with contributors Pedro Galeano, Alicia Nieto, and Eduardo Garcia-Portugues. See `NOTICE.md` and `upstream/DESCRIPTION` for provenance.

The translation focuses on numerical functional-data analysis kernels and uses array-oriented APIs rather than reproducing R S3 classes, formulas, plotting, or interactive behavior. Shared numerical infrastructure is reused through sibling FPM path dependencies on `rfortran-core` and `rfortran-linalg`; no BLAS, LAPACK, or dependency source is vendored here.

## Build

From this top-level package directory inside `Fortran-from-R-packages`:

```text
fpm build
fpm test
fpm run --example basic
```

The dependency layout expected by `fpm.toml` is:

```text
../rfortran-core
../rfortran-linalg
```

## Main numerical areas

The public `fdausc` module re-exports the package API. Implemented areas include kernel and integrated-kernel functions, Simpson integration and functional inner products, functional and multivariate distances, distance correlation, nonparametric smoother matrices and CV/GCV criteria, functional depths, functional PCA and PLS cores, conditional distribution/quantile/mode estimation, functional k-means and native k-NN/kernel classification, nonparametric/component regression cores, Gaussian and Ornstein-Uhlenbeck process simulation, accuracy/FDR helpers, and several functional-data test statistics.

## Translation coverage

Package status: **substantial**. Coverage is **179 of 179 (100.0%)** exported computational R functions. This fraction measures functions with a meaningful numerical mapping; it does **not** mean 100.0% R-interface compatibility. Formula parsing, S3/fdata object construction, plotting, smoothing-parameter estimation, resampling orchestration, and external-package dispatch remain different from R even where the corresponding numerical core is mapped.

The coverage denominator excludes 48 exported constructors/accessors/presentation helpers/R-only formula interfaces or functions that merely dispatch to an external package. `fpm.toml` records an empty untranslated list because all functions in this computational denominator have a mapping.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `AKer.cos` | `fdausc_kernels` | `aker_cos` | complete |
| `AKer.epa` | `fdausc_kernels` | `aker_epa` | complete |
| `AKer.norm` | `fdausc_kernels` | `aker_norm` | complete |
| `AKer.quar` | `fdausc_kernels` | `aker_quar` | complete |
| `AKer.tri` | `fdausc_kernels` | `aker_tri` | complete |
| `AKer.unif` | `fdausc_kernels` | `aker_unif` | complete |
| `Adot` | `fdausc_statistics` | `adot_matrix_vector` | complete |
| `CV.S` | `fdausc_smoothing` | `cv_score` | substantial |
| `FDR` | `fdausc_statistics` | `fdr_reject` | complete |
| `Ftest.statistic` | `fdausc_testing` | `flm_f_statistic` | complete |
| `GCCV.S` | `fdausc_smoothing` | `gccv_score` | substantial |
| `GCV.S` | `fdausc_smoothing` | `gcv_score` | substantial |
| `IKer.cos` | `fdausc_kernels` | `iker_cos` | complete |
| `IKer.epa` | `fdausc_kernels` | `iker_epa` | complete |
| `IKer.norm` | `fdausc_kernels` | `iker_norm` | complete |
| `IKer.quar` | `fdausc_kernels` | `iker_quar` | complete |
| `IKer.tri` | `fdausc_kernels` | `iker_tri` | complete |
| `IKer.unif` | `fdausc_kernels` | `iker_unif` | complete |
| `Ker.cos` | `fdausc_kernels` | `ker_cos` | complete |
| `Ker.epa` | `fdausc_kernels` | `ker_epa` | complete |
| `Ker.norm` | `fdausc_kernels` | `ker_norm` | complete |
| `Ker.quar` | `fdausc_kernels` | `ker_quar` | complete |
| `Ker.tri` | `fdausc_kernels` | `ker_tri` | complete |
| `Ker.unif` | `fdausc_kernels` | `ker_unif` | complete |
| `Kernel` | `fdausc_kernels` | `kernel_value` | complete |
| `Kernel.asymmetric` | `fdausc_kernels` | `asymmetric_kernel_value` | complete |
| `Kernel.integrate` | `fdausc_kernels` | `integrated_kernel_value` | complete |
| `LMDC.select` | `fdausc_extended` | `local_distance_correlation_select` | substantial |
| `MMD.test` | `fdausc_testing` | `mmd_rbf_statistic` | partial |
| `MMDA.test` | `fdausc_testing` | `mmd_rbf_statistic` | partial |
| `P.penalty` | `fdausc_smoothing` | `penalty_matrix` | partial |
| `PCvM.statistic` | `fdausc_statistics` | `pcvm_statistic_value` | complete |
| `S.KNN` | `fdausc_smoothing` | `smoothing_knn` | substantial |
| `S.LCR` | `fdausc_smoothing` | `smoothing_lcr` | substantial |
| `S.LLR` | `fdausc_smoothing` | `smoothing_llr` | substantial |
| `S.LPR` | `fdausc_smoothing` | `smoothing_lpr` | substantial |
| `S.NW` | `fdausc_smoothing` | `smoothing_nw` | substantial |
| `S.basis` | `fdausc_smoothing` | `basis_smoothing_matrix` | substantial |
| `Var.e` | `fdausc_smoothing` | `residual_variance` | substantial |
| `Var.y` | `fdausc_smoothing` | `fitted_variance` | substantial |
| `XYRP.test` | `fdausc_generalized` | `projected_two_sample_ks` | partial |
| `bcdcor.dist` | `fdausc_metrics` | `bias_corrected_distance_correlation` | substantial |
| `cat2meas` | `fdausc_statistics` | `confusion_accuracy, confusion_kappa, confusion_iou, confusion_recall, confusion_precision, confusion_f1` | partial |
| `classif.DD` | `fdausc_extended` | `depth_classify` | partial |
| `classif.depth` | `fdausc_extended` | `depth_classify` | partial |
| `classif.gkam` | `fdausc_generalized` | `multiclass_additive_classify` | partial |
| `classif.glm` | `fdausc_generalized` | `multiclass_glm_classify` | partial |
| `classif.gsam` | `fdausc_generalized` | `multiclass_penalized_glm_classify` | partial |
| `classif.gsam.vs` | `fdausc_generalized` | `forward_glm_select, multiclass_penalized_glm_classify` | partial |
| `classif.kernel` | `fdausc_models` | `kernel_classify` | substantial |
| `classif.kfold` | `fdausc_extended` | `select_classification_cv` | partial |
| `classif.knn` | `fdausc_models` | `knn_classify` | substantial |
| `classif.np` | `fdausc_models, fdausc_regression` | `kernel_classify, select_bandwidth_cv` | partial |
| `cond.F` | `fdausc_models` | `conditional_cdf` | substantial |
| `cond.mode` | `fdausc_models` | `conditional_mode` | substantial |
| `cond.quantile` | `fdausc_models` | `conditional_quantile` | substantial |
| `cov.test.fdata` | `fdausc_testing` | `covariance_hs_statistic` | partial |
| `dcor.dist` | `fdausc_metrics` | `distance_correlation` | substantial |
| `dcor.test` | `fdausc_metrics` | `distance_correlation_test` | substantial |
| `dcor.xy` | `fdausc_metrics` | `ordinary_distance_matrix, distance_correlation_test` | partial |
| `depth.FM` | `fdausc_depth` | `fraiman_muniz_depth` | substantial |
| `depth.FMp` | `fdausc_extended` | `multicomponent_mahalanobis_depth` | partial |
| `depth.FSD` | `fdausc_depth` | `functional_spatial_depth` | substantial |
| `depth.KFSD` | `fdausc_depth` | `kernel_functional_spatial_depth` | substantial |
| `depth.RP` | `fdausc_depth` | `projection_tukey_depth` | partial |
| `depth.RPD` | `fdausc_extended` | `rpd_depth_scores` | partial |
| `depth.RPp` | `fdausc_extended` | `rpd_depth_scores` | partial |
| `depth.RT` | `fdausc_depth` | `projection_tukey_depth` | partial |
| `depth.mode` | `fdausc_depth` | `depth_mode_scores` | substantial |
| `depth.modep` | `fdausc_extended` | `multicomponent_modal_depth` | partial |
| `dev.S` | `fdausc_generalized` | `deviance_smoother_score` | substantial |
| `dfv.statistic` | `fdausc_testing` | `dfv_statistics` | substantial |
| `dfv.test` | `fdausc_testing` | `dfv_statistics` | partial |
| `dis.cos.cor` | `fdausc_integration` | `cosine_proximity` | substantial |
| `fEqDistrib.test` | `fdausc_testing` | `distribution_equality_statistic` | partial |
| `fanova.RPm` | `fdausc_generalized` | `random_projection_anova` | partial |
| `fanova.hetero` | `fdausc_extended` | `hetero_anova_onefactor` | partial |
| `fanova.onefactor` | `fdausc_testing` | `fanova_mean_distance_statistic` | partial |
| `fdata.bootstrap` | `fdausc_extended` | `bootstrap_mean_replicates` | partial |
| `fdata.cen` | `fdausc_statistics` | `center_curves` | substantial |
| `fdata.deriv` | `fdausc_statistics` | `difference_derivative` | substantial |
| `fdata2basis` | `fdausc_extended` | `basis_project_grid` | substantial |
| `fdata2pc` | `fdausc_models` | `functional_pca` | substantial |
| `fdata2pls` | `fdausc_models` | `functional_pls1` | substantial |
| `flm.Ftest` | `fdausc_testing` | `flm_f_statistic` | partial |
| `flm.test` | `fdausc_statistics, fdausc_simulation` | `adot_matrix_vector, pcvm_statistic_value, wild_residuals` | partial |
| `fmean.test.fdata` | `fdausc_testing` | `mean_difference_statistic` | partial |
| `fregre.basis` | `fdausc_generalized` | `basis_scalar_regression` | partial |
| `fregre.basis.cv` | `fdausc_generalized, fdausc_extended` | `basis_scalar_regression, select_basis_gcv` | partial |
| `fregre.basis.fr` | `fdausc_generalized` | `basis_response_regression` | partial |
| `fregre.bootstrap` | `fdausc_generalized` | `bootstrap_linear_refits` | partial |
| `fregre.gkam` | `fdausc_generalized` | `additive_backfit` | partial |
| `fregre.glm` | `fdausc_generalized` | `glm_irls` | partial |
| `fregre.glm.vs` | `fdausc_generalized` | `forward_glm_select` | partial |
| `fregre.gls` | `fdausc_generalized` | `gls_fit` | partial |
| `fregre.gsam` | `fdausc_generalized` | `penalized_glm_irls` | partial |
| `fregre.gsam.vs` | `fdausc_generalized` | `forward_glm_select` | partial |
| `fregre.igls` | `fdausc_generalized` | `gls_fit` | partial |
| `fregre.lm` | `fdausc_regression` | `weighted_linear_fit` | partial |
| `fregre.np` | `fdausc_regression` | `nonparametric_regression` | substantial |
| `fregre.np.cv` | `fdausc_regression` | `select_bandwidth_cv` | partial |
| `fregre.pc` | `fdausc_models, fdausc_regression` | `functional_pca, component_regression` | partial |
| `fregre.pc.cv` | `fdausc_models, fdausc_regression` | `functional_pca, component_regression` | partial |
| `fregre.plm` | `fdausc_generalized` | `partially_linear_fit` | partial |
| `fregre.pls` | `fdausc_models, fdausc_regression` | `functional_pls1, component_regression` | partial |
| `fregre.pls.cv` | `fdausc_models, fdausc_regression` | `functional_pls1, component_regression` | partial |
| `func.mean` | `fdausc_statistics` | `functional_mean` | substantial |
| `func.med.FM` | `fdausc_depth` | `fraiman_muniz_depth, depth_summary` | partial |
| `func.med.RP` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.med.RPD` | `fdausc_extended, fdausc_depth` | `rpd_depth_scores, depth_summary` | partial |
| `func.med.RT` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.med.mode` | `fdausc_depth` | `depth_mode_scores, depth_summary` | partial |
| `func.trim.FM` | `fdausc_depth` | `fraiman_muniz_depth, depth_summary` | partial |
| `func.trim.RP` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.trim.RPD` | `fdausc_extended, fdausc_depth` | `rpd_depth_scores, depth_summary` | partial |
| `func.trim.RT` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.trim.mode` | `fdausc_depth` | `depth_mode_scores, depth_summary` | partial |
| `func.trimvar.FM` | `fdausc_depth` | `fraiman_muniz_depth, depth_summary` | partial |
| `func.trimvar.RP` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.trimvar.RPD` | `fdausc_extended, fdausc_depth` | `rpd_depth_scores, depth_summary` | partial |
| `func.trimvar.RT` | `fdausc_depth` | `projection_tukey_depth, depth_summary` | partial |
| `func.trimvar.mode` | `fdausc_depth` | `depth_mode_scores, depth_summary` | partial |
| `func.var` | `fdausc_statistics` | `functional_variance` | substantial |
| `gridfdata` | `fdausc_simulation` | `grid_basis_combinations` | substantial |
| `h.default` | `fdausc_extended` | `default_bandwidth_grid` | substantial |
| `influence_quan` | `fdausc_generalized` | `influence_quantile_summary` | partial |
| `inprod.fdata` | `fdausc_integration` | `inprod_curves` | substantial |
| `int.simpson` | `fdausc_integration` | `integrate_curve` | substantial |
| `int.simpson2` | `fdausc_integration` | `integrate_curve` | substantial |
| `kmeans.center.ini` | `fdausc_models` | `functional_kmeans` | partial |
| `kmeans.fd` | `fdausc_models` | `functional_kmeans` | substantial |
| `ldata.cen` | `fdausc_statistics` | `center_curves` | partial |
| `mdepth.FM` | `fdausc_depth` | `tukey_coordinate_depth` | partial |
| `mdepth.FSD` | `fdausc_depth` | `multivariate_spatial_depth` | substantial |
| `mdepth.HS` | `fdausc_depth` | `projection_tukey_depth` | partial |
| `mdepth.KFSD` | `fdausc_depth` | `multivariate_kernel_spatial_depth` | substantial |
| `mdepth.LD` | `fdausc_extended` | `likelihood_depth` | substantial |
| `mdepth.MhD` | `fdausc_depth` | `mahalanobis_depth` | substantial |
| `mdepth.RP` | `fdausc_depth` | `projection_tukey_depth` | partial |
| `mdepth.SD` | `fdausc_extended` | `bivariate_simplicial_depth` | partial |
| `mdepth.TD` | `fdausc_depth` | `tukey_coordinate_depth` | substantial |
| `metric.DTW` | `fdausc_metrics` | `dtw_distance_matrix` | substantial |
| `metric.TWED` | `fdausc_metrics` | `twed_distance_matrix` | substantial |
| `metric.WDTW` | `fdausc_metrics` | `wdtw_distance_matrix` | substantial |
| `metric.dist` | `fdausc_metrics` | `ordinary_distance_matrix` | substantial |
| `metric.hausdorff` | `fdausc_metrics` | `hausdorff_distance_matrix` | substantial |
| `metric.kl` | `fdausc_metrics` | `kl_distance_matrix` | substantial |
| `metric.ldata` | `fdausc_metrics` | `combine_distance_matrices` | partial |
| `metric.lp` | `fdausc_metrics` | `lp_distance_matrix` | substantial |
| `metric.mfdata` | `fdausc_metrics` | `combine_distance_matrices` | partial |
| `mfdata.cen` | `fdausc_statistics` | `center_curves` | partial |
| `norm.fdata` | `fdausc_integration` | `norm_curves` | substantial |
| `optim.basis` | `fdausc_extended` | `select_basis_gcv` | partial |
| `optim.np` | `fdausc_regression` | `select_bandwidth_cv` | partial |
| `outliers.depth.pond` | `fdausc_extended` | `depth_outlier_flags` | partial |
| `outliers.depth.trim` | `fdausc_extended` | `depth_outlier_flags` | partial |
| `outliers.lrt` | `fdausc_extended` | `lrt_outlier_statistic` | partial |
| `outliers.thres.lrt` | `fdausc_extended` | `lrt_outlier_statistic` | partial |
| `pred.MAE` | `fdausc_statistics` | `pred_mae` | complete |
| `pred.MSE` | `fdausc_statistics` | `pred_mse` | complete |
| `pred.RMSE` | `fdausc_statistics` | `pred_rmse` | complete |
| `pred2meas` | `fdausc_statistics` | `prediction_measures` | complete |
| `pvalue.FDR` | `fdausc_statistics` | `fdr_adjusted_pvalue` | complete |
| `r.ou` | `fdausc_simulation` | `simulate_ou_process` | substantial |
| `rcombfdata` | `fdausc_simulation` | `random_basis_combinations` | substantial |
| `rdir.pc` | `fdausc_models, fdausc_simulation` | `functional_pca, random_basis_combinations` | partial |
| `rp.flm.statistic` | `fdausc_statistics` | `rp_projection_statistics` | complete |
| `rp.flm.test` | `fdausc_statistics, fdausc_simulation` | `rp_projection_statistics, wild_residuals` | partial |
| `rproc2fdata` | `fdausc_simulation` | `process_covariance, simulate_gaussian_process` | substantial |
| `rwild` | `fdausc_simulation` | `wild_residuals` | substantial |
| `semimetric.basis` | `fdausc_extended` | `basis_semimetric` | partial |
| `semimetric.deriv` | `fdausc_statistics, fdausc_models` | `difference_derivative, derivative_distance` | substantial |
| `semimetric.fourier` | `fdausc_extended` | `fourier_semimetric` | partial |
| `semimetric.hshift` | `fdausc_extended` | `horizontal_shift_distance_matrix` | substantial |
| `semimetric.mplsr` | `fdausc_models` | `functional_pls1, pca_score_distance` | partial |
| `semimetric.pca` | `fdausc_models` | `functional_pca, pca_score_distance` | partial |
| `tab2meas` | `fdausc_statistics` | `confusion_accuracy, confusion_kappa, confusion_iou, confusion_iou_by_class, confusion_recall, confusion_precision, confusion_f1` | substantial |
| `trace.matrix` | `fdausc_smoothing` | `matrix_trace` | substantial |
| `weights4class` | `fdausc_statistics` | `class_weights` | substantial |

## Compatibility notes

`real(dp)` is used consistently, with `dp` imported from `rfortran-core` (`r_kinds`). The translation intentionally does not seek bit-for-bit RNG parity with R. Procedures that need linear solves, eigendecompositions, or SVD use the public `rfortran-linalg` API rather than system BLAS/LAPACK links. Tests are deterministic and avoid depending on random draws.

See `docs/API_COVERAGE.md` for the omitted-family summary and `docs/VALIDATION.md` for validation details.
