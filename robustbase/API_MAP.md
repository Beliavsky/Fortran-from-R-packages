# API map

This map distinguishes definition-preserving translations, tested numerical reimplementations, and excluded R infrastructure.

## Scales, scores, and descriptive utilities

| robustbase family | Fortran API | Status |
|---|---|---|
| `wgt.himedian` | `weighted_high_median` | Strict-half weighted high median; tested |
| `Qn`, `Sn` | `qn_scale`, `sn_scale` | Exact pairwise definitions and default corrections; tested |
| `s_mad`, `s_IQR` | `mad_scale`, `iqr_scale` | Tested |
| `huberM`, `huberize`, `tauHuber` | `huber_location`, `huberize_vector`, `tau_huber` | Tested |
| `scaleTau2` | `scale_tau2` | Tested numerical implementation |
| Huber/Hampel/Tukey psi families | matching `*_psi`, `*_rho`, `*_weight` routines | Tested |
| Welsh, optimal, GGW, LQQ psi families | matching `welsh_*`, `optimal_*`, `ggw_*`, `lqq_*` routines | Score/loss/weight derivatives tested; GGW rho uses numerical integration |
| `mc`, `lmc`, `rmc` | `medcouple`, `left_medcouple`, `right_medcouple` | Exact naive kernel median; tested |
| `adjboxStats` | `adjusted_boxplot_stats` | Tested |
| `rowMedians`, `colMedians` | `row_medians`, `column_medians` | Tested |
| `rankMM`, `classPC`, `.signflip` | `rank_mm`, `classical_pca` | Tested |
| `fullRank` | `full_rank_matrix`, `independent_columns` | Tested |

## Covariance and outlyingness

| robustbase family | Fortran API | Status |
|---|---|---|
| comedian / `COM` | `comedian`, `cov_comedian` | Tested |
| `covComed` | `cov_comed` | Iterative numerical reimplementation; tested |
| `covGK` | `cov_gk` | Tested |
| `covOGK` | `cov_ogk` | Fixed Qn/GK implementation; tested |
| `covMcd` basic | `cov_mcd` | Random-start C-step estimator; tested |
| `detMCD`, `r6pack` | `cov_detmcd` | Six starts, C-steps, reweighting, exact-fit hyperplane output; tested |
| FAST-MCD partition/refine path | `fast_mcd_partitioned` | Tested numerical reimplementation |
| `.MCDcons`, `.MCDcnp2`, `.MCDcnp2.rew` | `mcd_consistency_factor`, `mcd_finite_sample_factor`, `mcd_reweighted_finite_sample_factor` | Analytical formulas; tested |
| `robMD`, `mahalanobisD` | `robust_mahalanobis` | Tested |
| projection outlyingness | `adjusted_outlyingness` | Tested |
| `adjOutlyingness` continuous path | `adjusted_outlyingness_full` | Hyperplane-direction generation and skew-adjusted distances; tested |
| `tolEllipsePlot` numerical core | `tolerance_ellipse_points` | Coordinates only; tested |

## Linear regression

| robustbase family | Fortran API | Status |
|---|---|---|
| `ltsReg` basic | `lts_regression` | Random-start C-step estimator; tested |
| `ltsReg` advanced | `fast_lts_regression` | Exact/best/deterministic/random starts and reweighting; tested |
| partitioned FAST-LTS | `fast_lts_partitioned` | Partition candidates plus full-data refinement; tested |
| `lmrob.lar` | `lmrob_lar_fit` | Smoothed LAD IRLS numerical implementation; tested |
| `lmrob.S` | `lmrob_s_fit` | S scale, nonsingular/exact/simple/best subsampling; tested |
| `lmrob.fit.MM`, `lmrob..M..fit`, `lmrob..D..fit` | `lmrob_fit` with `S`, `SM`/`MM`, or `SMDM` | Tested |
| basic `lmrob` analogue | `mm_regression` | LTS plus robust IRLS; tested |
| `predict.lmrob` numerical core | `robust_linear_predict` | Confidence/prediction intervals; tested |
| `anova.lmrob` numerical cores | `robust_wald_test`, `robust_deviance_test` | Tested |
| `outlierStats`, robust fit measure | `robust_outlier_stats`, `robust_r_squared` | Tested |

## GLM and nonlinear regression

| robustbase family | Fortran API | Status |
|---|---|---|
| basic robust GLM | `robust_glm_fit` | Huberized binomial/Poisson IRLS; tested |
| `BYlogreg` | `by_logistic_fit`, `by_phi*` | Objective, derivatives, robust updates, covariance; tested |
| `glmrob` Mqle | `glmrob_mqle_fit` | Binomial/Poisson bias-corrected robust estimating equations; tested |
| `glmrob` MT | `glmrob_mt_fit` | Binomial/grouped-binomial transformed estimator; tested |
| basic `nlrob` analogue | `robust_nls_fit` | Finite-difference robust IRLS; tested |
| `nlrob` MM/tau/CM/MTL | `nlrob_mm_fit`, `nlrob_tau_fit`, `nlrob_cm_fit`, `nlrob_mtl_fit` | Bounded numerical reimplementations; tested |

## Supporting APIs

- LAPACK-backed least squares, eigendecomposition, SVD rank, symmetric inversion, and covariance
- Normal and chi-square CDF/quantile utilities
- Two-column CSV reader
- Result types carrying estimates, residuals, weights, covariance matrices, standard errors, objectives, iteration counts, and convergence flags

## Excluded compatibility surfaces

- R formulas, model frames, factors/contrasts, offsets, family objects, and missing-value policies
- S3/S4 classes and methods, summaries, printing, plotting, and package data
- Categorical adjusted-outlyingness logic
- Exact historical finite-sample lookup-table overrides
- Efficiency/breakdown-point calibration solvers for GGW and LQQ tuning constants
- Exact iteration and random-stream equivalence with the legacy optimized kernels
