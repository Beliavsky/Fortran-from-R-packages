# Selected API mapping

| R / C source | Fortran API |
|---|---|
| `clo` | `closure`, `closure_rows` |
| acomp `+`, `*` geometry | `perturb`, `power_comp` |
| `clr`, `clrInv` | `clr`, `clr_inv` |
| `ilr`, `ilrInv`, `ilrBase` | `ilr`, `ilr_inv`, `ilr_base` |
| `alr`, `alrInv` | `alr`, `alr_inv` |
| `apt`, `aptInv` | `apt`, `apt_inv` |
| `cpt`, `cptInv` | `cpt`, `cpt_inv` |
| `ipt`, `iptInv` | `ipt`, `ipt_inv` |
| `ilt`, `iltInv` | `ilt`, `ilt_inv` |
| `variation` | `variation_matrix` |
| `clrvar2ilr`, `ilrvar2clr` | `clrvar_to_ilr`, `ilrvar_to_clr` |
| `clrvar2variation`, inverse | `clrvar_to_variation`, `variation_to_clrvar` |
| `mean.acomp` | `acomp_mean` |
| `fitDirichlet` | `fit_dirichlet` |
| `ddirichlet.acomp` | `dirichlet_pdf`, `dirichlet_logpdf` |
| Dirichlet RNG | `rdirichlet` |
| `dnorm.acomp` | `logistic_normal_pdf`, `logistic_normal_logpdf` |
| `rnorm.acomp` | `rlogistic_normal` |
| `AitchisonDistributionIntegrals` | `aitchison_integrals` |
| `dAitchison` | `aitchison_pdf`, `aitchison_logpdf` |
| `rAitchison` default path | `raitchison` |
| classical `var.acomp` | `compositional_covariance` |
| robust `var.acomp` / `covMcd` path | `robust_compositional_covariance` |
| `princomp.acomp` / robust PCA core | `compositional_pca` |
| `acompNormalLocation.test` | `acomp_normal_location_one_sample`, two-sample routine |
| compositional `lm` numerical core | `compositional_lm_fit`, `compositional_lm_predict` |
| `gsi.PrinBal(..., "PBhclust")` | `principal_balance_hclust` |
| `gsi.PrinBal(..., "PBangprox")` | `principal_balance_angprox` |
| `gsi.PrinBal(..., "PBmaxvar")` | `principal_balance_maxvar` |
| `gsiCGSvariogram` complete-data core | `logratio_variogram` |
| `gsiCGSvg2lrvg` | `vg_to_lrvg` |
| variogram model functions | `vgram_*` |
| ordinary complete-data kriging core | `compositional_ordinary_kriging` |
| `gsiCGSkriging` predictor | `compositional_general_kriging` |
| `gsiCImpAcompGetTypes` / `GetIdx` | `missing_pattern_indices` |
| `gsiCImpAcompNewImputationVariance` | `conditional_alr_moments` |
| `gsiCImpAcompCompleteAlr` + BDL expectation core | `impute_acomp_conditional` |
| `gsiCImpAcompFitWithProjection` | `fit_acomp_projection` |
| completed iterative imputation extension | `fit_acomp_em` |
| `rmultinom.ccomp` | `rmultinom_composition` |
| Poisson count simulation | `rpois_composition` |
| `gsiDensityCheck` | `kernel_similarity_statistic`, permutations |
| `gsiKSPoisson` | `poisson_ks_statistic` |
| `gsiKSsortedUniforms` intended kernel | `sorted_uniforms` |
| `gsiKSPoissonSample` | `poisson_ks_sample` |
| `gsiCImpAcompClrExpectation` / cache setup | `build_imputation_cache`, `acomp_clr_expectation` |
| empirical/max Mahalanobis simulation | `r_empirical_mahalanobis`, `r_max_mahalanobis` |
| empirical/max Mahalanobis calibration | `q_*_mahalanobis`, `p_*_mahalanobis`, `is_mahalanobis_outlier` |
| `OutlierClassifier1` numerical core | `outlier_classifier_best` |
| `gsi.AcompGOFEtest` energy wrapper | `acomp_energy_test` |
| `acompNormalGOF.test` energy wrapper | `acomp_normal_energy_test` |
| tensorA named/indexed operations | `compositions_tensor` / compiled `tensora` modules |
