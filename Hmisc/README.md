# Hmisc - modern Fortran computational translation

This package translates a focused numerical subset of the R package **Hmisc 5.3-0** to modern free-form Fortran with FPM layout.  It preserves upstream provenance and GPL licensing while omitting plotting, interactive interfaces, R classes, formula processing, HTML/LaTeX/Typst output, and data-import machinery.

No external BLAS/LAPACK installation is required for the current subset.

## Build

```text
fpm build
fpm test
fpm run --example basic_stats
```

The implementation uses a single `dp` real kind based on `real64`, explicit interfaces, and `ieee_is_nan` for missing numeric values represented as IEEE NaNs.

## Translation coverage

Package status: **substantial**.

Package-wide computational coverage: **136 of 142 (95.8%)** exported computational R functions are mapped. This fraction counts mapped computational functions, not complete R interface compatibility. Plotting, printing, formatting, data-download, interactive-display, presentation-only exports, and internal/package-dispatch helpers are excluded from the denominator. The 6 currently untranslated computational exports are listed exactly in `fpm.toml`. The audited denominator excludes `ordGridFun` (graphics dispatch), `nFm` (numeric-string formatting), and `var.inner` (an internal formula-parsing helper) under the stated coverage rules.

The remaining six computational exports are intentionally deferred because their essential behavior is inseparable from large R fitting, imputation, transformation, clustering, or callback frameworks; they are not represented by placeholder mappings.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `GiniMd` | `hmisc` | `gini_md` | complete |
| `trap.rule` | `hmisc` | `trap_rule` | complete |
| `samplesize.bin` | `hmisc` | `samplesize_bin` | substantial |
| `wtd.mean` | `hmisc` | `weighted_mean` | substantial |
| `wtd.var` | `hmisc` | `weighted_variance` | substantial |
| `wtd.quantile` | `hmisc` | `weighted_quantile` | partial |
| `wtd.rank` | `hmisc` | `weighted_rank` | substantial |
| `wtd.table` | `hmisc` | `weighted_table` | partial |
| `somers2` | `hmisc` | `somers2` | substantial |
| `rcorr` | `hmisc` | `rcorr` | substantial |
| `hoeffd` | `hmisc` | `hoeffd` | substantial |
| `rcorr.cens` | `hmisc` | `rcorr_cens` | substantial |
| `rcorrp.cens` | `hmisc` | `rcorrp_cens` | substantial |
| `cutGn` | `hmisc` | `cut_gn` | partial |
| `whichClosest` | `hmisc` | `which_closest` | complete |
| `whichClosePW` | `hmisc` | `which_close_pw` | substantial |
| `pMedian` | `hmisc` | `pseudomedian` | partial |
| `dualSD` | `hmisc` | `dual_sd` | substantial |
| `bpower` | `hmisc` | `bpower` | substantial |
| `bsamsize` | `hmisc` | `bsamsize` | complete |
| `ballocation` | `hmisc` | `ballocation` | substantial |
| `Weibull2` | `hmisc` | `weibull2_fit`, `weibull_survival` | complete |
| `Gompertz2` | `hmisc` | `gompertz2_fit`, `gompertz_survival` | complete |
| `Lognorm2` | `hmisc` | `lognorm2_fit`, `lognorm_survival` | complete |
| `improveProb` | `hmisc` | `improve_prob` | partial |
| `rcspline.eval` | `hmisc` | `rcspline_eval` | partial |
| `popower` | `hmisc` | `popower` | complete |
| `posamsize` | `hmisc` | `posamsize` | complete |
| `pomodm` | `hmisc` | `pomodm` | partial |
| `groupn` | `hmisc` | `groupn` | substantial |
| `spearman` | `hmisc` | `spearman_rho` | substantial |
| `stepfun.eval` | `hmisc` | `stepfun_eval` | partial |
| `xySortNoDupNoNA` | `hmisc` | `xy_sort_no_dup_no_na` | substantial |
| `nCoincident` | `hmisc` | `n_coincident` | substantial |
| `matxv` | `hmisc` | `matxv` | partial |
| `cpower` | `hmisc` | `cpower` | substantial |
| `ciapower` | `hmisc` | `ciapower` | substantial |
| `wtd.Ecdf` | `hmisc` | `wtd_ecdf` | substantial |
| `smean.sd` | `hmisc` | `smean_sd` | complete |
| `smean.sdl` | `hmisc` | `smean_sdl` | complete |
| `smedian.hilow` | `hmisc` | `smedian_hilow` | substantial |
| `approxExtrap` | `hmisc` | `approx_extrap` | partial |
| `james.stein` | `hmisc` | `james_stein` | substantial |
| `ftupwr` | `hmisc` | `ftupwr` | complete |
| `ftuss` | `hmisc` | `ftuss` | complete |
| `ecdfSteps` | `hmisc` | `ecdf_steps` | substantial |
| `Lag` | `hmisc` | `lag_numeric` | partial |
| `nomiss` | `hmisc` | `nomiss` | substantial |
| `fillin` | `hmisc` | `fillin` | substantial |
| `cumcategory` | `hmisc` | `cumcategory` | substantial |
| `deff` | `hmisc` | `deff` | substantial |
| `R2Measures` | `hmisc` | `r2_measures` | substantial |
| `binconf` | `hmisc` | `binconf_wilson`, `binconf_asymptotic` | partial |
| `invertTabulated` | `hmisc` | `invert_tabulated_linear` | partial |
| `score.binary` | `hmisc` | `score_binary_max`, `score_binary_sum` | partial |
| `mhgr` | `hmisc` | `mhgr` | substantial |
| `lrcum` | `hmisc` | `lrcum` | substantial |
| `num.denom.setup` | `hmisc` | `num_denom_setup` | complete |
| `rMultinom` | `hmisc` | `rmultinom` | substantial |
| `qrxcenter` | `hmisc` | `qrxcenter` | substantial |
| `solvet` | `hmisc` | `solvet`, `solvet_inverse` | substantial |
| `pc1` | `hmisc` | `pc1` | substantial |
| `smean.cl.normal` | `hmisc` | `smean_cl_normal` | substantial |
| `whichClosek` | `hmisc` | `which_close_k` | substantial |
| `lm.fit.qr.bare` | `hmisc` | `lm_fit_qr_bare` | substantial |
| `abs.error.pred` | `hmisc` | `abs_error_pred` | substantial |
| `rcorrcens` | `hmisc` | `rcorrcens_summary` | partial |
| `logrank` | `hmisc` | `logrank` | substantial |
| `spearman2` | `hmisc` | `spearman2_numeric` | partial |
| `equalBins` | `hmisc` | `equal_bins` | substantial |
| `jitter2` | `hmisc` | `jitter2_numeric` | partial |
| `hdquantile` | `hmisc` | `hdquantile` | partial |
| `catTestchisq` | `hmisc` | `cat_test_chisq` | substantial |
| `conTestkw` | `hmisc` | `con_test_kw` | substantial |
| `t.test.cluster` | `hmisc` | `t_test_cluster` | substantial |
| `km.quick` | `hmisc` | `km_quick` | substantial |
| `bystats` | `hmisc` | `bystats_mean` | partial |
| `bystats2` | `hmisc` | `bystats2_mean` | partial |
| `all.digits` | `hmisc` | `all_digits` | complete |
| `all.is.numeric` | `hmisc` | `all_is_numeric` | partial |
| `find.matches` | `hmisc` | `find_matches` | substantial |
| `seqFreq` | `hmisc` | `seq_freq` | substantial |
| `pairUpDiff` | `hmisc` | `pair_up_diff` | substantial |
| `smean.cl.boot` | `hmisc` | `smean_cl_boot` | substantial |
| `smearingEst` | `hmisc` | `smearing_est_tabulated` | partial |
| `bpower.sim` | `hmisc` | `bpower_sim_counts` | partial |
| `spower` | `hmisc` | `spower_simulated` | partial |
| `inverseFunction` | `hmisc` | `inverse_function_all` | partial |
| `bezier` | `hmisc` | `bezier_curve` | substantial |
| `cut2` | `hmisc` | `cut2_explicit` | partial |
| `simPOcuts` | `hmisc` | `sim_po_cuts` | substantial |
| `rcsplineFunction` | `hmisc` | `rcspline_function` | substantial |
| `rcspline.restate` | `hmisc` | `rcspline_restate_coefficients` | partial |
| `spearman.test` | `hmisc` | `spearman_test` | substantial |
| `chiSquare` | `hmisc` | `chi_square_codes` | partial |
| `combine.levels` | `hmisc` | `combine_levels_codes` | substantial |
| `princmp` | `hmisc` | `princmp_regular` | partial |
| `bootkm` | `hmisc` | `bootkm_resampled` | substantial |
| `matchCases` | `hmisc` | `match_cases_core` | substantial |
| `largest.empty` | `hmisc` | `largest_empty_rexhaustive` | partial |
| `gbayes` | `hmisc` | `gbayes_update` | substantial |
| `gbayesMixPredNoData` | `hmisc` | `gbayes_mix_pred_density`, `gbayes_mix_pred_cdf` | substantial |
| `gbayesMixPost` | `hmisc` | `gbayes_mix_post_density`, `gbayes_mix_post_cdf`, `gbayes_mix_post_mean` | substantial |
| `gbayes1PowerNP` | `hmisc` | `gbayes1_power_np` | complete |
| `gbayesMixPowerNP` | `hmisc` | `gbayes_mix_power_np` | partial |

| `gbayes2` | `hmisc` | `gbayes2_tabulated` | partial |
| `propsPO` | `hmisc` | `props_po_counts` | partial |
| `propsTrans` | `hmisc` | `props_trans_counts` | substantial |
| `ynbind` | `hmisc` | `ynbind_numeric` | partial |
| `reformM` | `hmisc` | `reformm_missing_order` | partial |
| `mApply` | `hmisc` | `mapply_group_mean` | partial |
| `soprobMarkovOrd` | `hmisc` | `soprob_markov_ord` | substantial |
| `ordGroupBoot` | `hmisc` | `ord_group_boot_mean` | substantial |
| `clowess` | `hmisc` | `clowess_smooth` | substantial |
| `curveSmooth` | `hmisc` | `curve_smooth_observed` | partial |
| `wtd.loess.noiter` | `hmisc` | `wtd_loess_noiter` | substantial |
| `simMarkovOrd` | `hmisc` | `sim_markov_ord` | substantial |
| `impute` | `hmisc` | `impute_median`, `impute_constant` | partial |
| `movStats` | `hmisc` | `mov_stats_n_raw` | partial |
| `summarize` | `hmisc` | `summarize_mean_codes` | partial |
| `aregTran` | `hmisc` | `areg_tran_numeric` | substantial |
| `areg` | `hmisc` | `areg_linear_fit` | partial |
| `dataframeReduce` | `hmisc` | `dataframe_reduce_numeric` | partial |
| `dataRep` | `hmisc` | `data_rep_numeric` | partial |
| `varclus` | `hmisc` | `varclus_numeric` | substantial |
| `redun` | `hmisc` | `redun_numeric` | partial |
| `gbayesSeqSim` | `hmisc` | `gbayes_seq` | substantial |
| `areg.boot` | `hmisc` | `areg_boot_linear` | partial |
| `fit.mult.impute` | `hmisc` | `fit_mult_impute_combine` | substantial |
| `impute.transcan` | `hmisc` | `impute_transcan_numeric` | partial |
| `estSeqSim` | `hmisc` | `est_seq_sim_mean_difference` | partial |
| `rm.boot` | `hmisc` | `rm_boot_linear` | partial |
| `ordTestpo` | `hmisc` | `ord_test_po_numeric` | substantial |
| `intMarkovOrd` | `hmisc` | `int_markov_ord_intercepts` | partial |
| `simRegOrd` | `hmisc` | `sim_reg_ord_supplied` | partial |
| `soprobMarkovOrdm` | `hmisc` | `soprob_markov_ordm_core` | substantial |

See `API_COVERAGE.md` and the `[[extra.translation.function]]` records in `fpm.toml` for compatibility notes.

## Validation

The deterministic test program covers weighted statistics, correlations, Hoeffding D, censored concordance, grouping, nearest-value lookup, pseudomedian, Gini mean difference, trapezoidal integration, binary sample-size/power/allocation calculations, dual standard deviations, two-point survival-distribution fitting, continuous NRI components, restricted cubic spline basis evaluation, proportional-odds power/sample-size calculations, rank/step helpers, coincident-point counting, and matrix-vector prediction helpers, survival-study power, weighted ECDFs, descriptive interval summaries, linear interpolation/extrapolation, James-Stein group shrinkage, Fleiss-Tytun-Ury binary power/sample-size approximations, ECDF coordinates, numeric lagging, missing-value filtering/filling, cumulative-category encoding, cluster design effects, generalized R-squared measures, binomial confidence intervals, tabulated-function inversion, binary-item scoring, Mantel-Haenszel risk ratios, cumulative diagnostic likelihood ratios, numerator/denominator expansion, multinomial sampling with supplied uniforms, centered QR transforms, QR least-squares/inversion, first-principal-component scoring, normal-theory confidence summaries, and deterministic k-nearest selection, Kaplan-Meier survival estimation, grouped default means, character numeric/digit tests, numeric tolerance matching, hierarchical sequential-condition assignment, and paired-difference confidence calculations, bootstrap percentile confidence limits for means, tabulated-inverse smearing estimates, deterministic binary-power simulation aggregation, supplied-time survival power simulation, and piecewise multi-branch function inversion, Bezier evaluation, explicit-cut categorization, deterministic proportional-odds cut simulation, and restricted-cubic-spline function/restate numerical cores, rank-regression Spearman tests, categorical pooling/chi-square testing, regular numeric PCA, and Gaussian conjugate/mixture Bayesian calculations and mixture-posterior power root solving, robust LOWESS, weighted local-polynomial smoothing, bootstrap-safe ordinal grouping, per-curve observed-point smoothing, deterministic ordinal Markov simulation, numeric imputation, raw fixed-sample moving summaries, multi-key grouped means, areg design transformations and linear fitting, numeric data-frame reduction/representativeness, variable clustering/redundancy elimination, sequential Gaussian posterior assertions, and fitted-Markov transition occupancy propagation.  The code is also tested with gfortran bounds/runtime checking enabled.

## License and provenance

See `NOTICE.md` and `LICENSES/GPL-2.0-or-later.txt`.
