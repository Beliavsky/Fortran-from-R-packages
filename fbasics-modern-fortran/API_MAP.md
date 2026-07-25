# API map

This map describes the computational surface in version 0.3.0. Plain-array Fortran procedures and derived result types replace R formula, S3/S4, and time-index interfaces.

| fBasics area or routine | Modern Fortran procedure(s) | Status |
|---|---|---|
| `hilbert`, `pascal` | `hilbert_matrix`, `pascal_matrix` | Implemented and tested |
| `inv`, `rk`, matrix norms, `tr` | `matrix_inverse`, `matrix_rank`, `matrix_norm`, `matrix_trace` | Implemented and tested |
| `kron`, `triang`, `vec`, `vech` | `kronecker_product`, triangle procedures, `vec_matrix`, `vech_matrix` | Implemented and tested |
| `isPositiveDefinite`, `makePositiveDefinite` | `is_positive_definite`, `make_positive_definite` | Implemented and tested |
| `gridVector`, `tslag`, `pdl` | `grid_vector`, `lag_matrix`, `polynomial_distributed_lags` | Implemented and tested |
| `basicStats`, row statistics | `basic_stats`, `row_*` procedures | Implemented and tested |
| Robust sample moments and `sampleLmoments` | `sample_iqr`, robust skew/kurtosis, `sample_lmoments` | Implemented and tested |
| LCG and basic random generators | `set_lcg_seed`, `runif_lcg`, `rnorm_lcg`, `rt_lcg`, Gamma/chi-square/inverse-Gaussian procedures | Implemented and tested |
| Normal distribution and `nFit` | `dnorm_fs`, `pnorm_fs`, `qnorm_fs`, `rnorm_fs`, `fit_normal` | Implemented and tested |
| Student distribution and `tFit` | `dt_fs`, `pt_fs`, `qt_fs`, `rt_fs`, `fit_student` | Implemented and tested |
| `stableFit` and stable-law calculations | `dstable_s1`, `pstable_s1`, `qstable_s1`, `rstable_s1`, `fit_stable_ecf`, `fit_stable_mle` | Numerical S1 implementation, tested |
| NIG distribution and `nigFit` | NIG density/CDF/quantile/RNG/moments and `fit_nig` | Implemented and tested |
| GH/HYP/GHT distribution calculations | GH-family density/CDF/quantile/RNG/moment procedures | Implemented and tested |
| GH/HYP/GHT fitting wrappers | `fit_gh`, `fit_hyp`, `fit_ght` | Numerical analogues, tested |
| SGH/SNIG/SGHT wrappers | standardized procedures plus `fit_sgh`, `fit_snig`, `fit_sght` | Implemented and tested |
| GH-family robust moments | `gh_robust_moments`, `hyp_robust_moments`, `ght_robust_moments`, `nig_robust_moments`, `sgh_robust_moments`, `snig_robust_moments` | Implemented and tested |
| Ramberg-Schmeiser GLD | `dgld_rs`, `pgld_rs`, `qgld_rs`, `rgld_rs`, `gld_mode`, `fit_gld_quantiles` | Implemented and tested |
| FMKL GLD | `dgld_fmkl`, `pgld_fmkl`, `qgld_fmkl`, `rgld_fmkl` | Implemented and tested |
| Five-parameter GLD | `dgld_fm5`, `pgld_fm5`, `qgld_fm5`, `rgld_fm5` | Implemented and tested |
| GLD fitting modes | `fit_gld_extended` with `rob`, `mle`, `mps`, `gof`, `hist` | Numerical analogues, tested |
| `ssdFit`, `dssd`, `pssd`, `qssd`, `rssd` | `fit_spline_density`, `dssd`, `pssd`, `qssd`, `rssd` | Penalized B-spline analogue, tested |
| `.gmm` | `fit_gmm` identity/two-step/iterated/CUE | Implemented and tested |
| `.gel` | `fit_gel` EL/ET/CUE/ETEL | Implemented and tested |
| `.HAC`, `.kweights`, bandwidths | `moment_covariance`, five kernels, `newey_west_bandwidth`, `andrews_bandwidth`, `andrews_bandwidth_value` | AR(1) and ARMA(1,1) plug-in paths implemented and tested |
| GMM prewhitening | `prewhiten_var`; `fit_gmm(..., prewhite_order=...)` | VAR(p) residual prewhitening and long-run recoloring implemented; orders 1 and 2 tested |
| ARMA bandwidth approximation | `fit_arma11_css`; `bandwidth_method="andrews_arma11"` | Conditional-sum-of-squares numerical analogue, tested |
| GMM J test and `.lintest` | `gmm_result` J/bandwidth/prewhitening fields, `linear_restriction_test` | Implemented and tested |
| `linearInterp`, `linearInterpp` | `linear_interp`, `bilinear_interp` | Implemented and tested |
| `akimaInterp`, `akimaInterpp` | `triangulated_interp`, `local_plane_interp` | Tested linear numerical alternatives; not exact Akima cubic |
| `krigeInterp` | `estimate_kriging_model`, `ordinary_kriging` | Ordinary-kriging numerical implementation, tested |
| Correlation tests | Pearson, Spearman, Kendall procedures | Implemented and tested |
| Location and variance tests | `location_t_test`, `variance_f_test` | Implemented and tested |
| Jarque-Bera LM/ALM area | `jarque_bera_test`, `adjusted_jarque_bera_test` | Formula-based calculations tested; historical lookup tables omitted |
| Shapiro-Wilk/Francia and other normality tests | Corresponding public test procedures | Implemented and tested |
| Ansari-Bradley, Mood, Bartlett, Fligner-Killeen | Corresponding public procedures | Implemented and tested |
| Wilcoxon and Kruskal-Wallis | Corresponding public procedures | Implemented and tested |
| `ks2Test` and one-sample normal tests | KS/CvM/AD/Lilliefors/Pearson procedures | Implemented with documented approximations, tested |
| `maxddStats`, `dmaxdd`, `pmaxdd`, `rmaxdd` | Drawdown procedures | Implemented and tested |
| `tsHessian` | `numerical_hessian` | Numerical analogue, exercised by fits |
| Plot/GUI/color/S3/S4/formula/time-index code | None | Excluded as infrastructure |
