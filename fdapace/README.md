# fdapace: modern Fortran translation

This directory contains a modern free-form Fortran translation of computational routines from the R package **fdapace 0.6.0** (CRAN publication date 2024-07-03). The original package provides functional-data analysis and empirical-dynamics methods. This translation focuses on reusable numerical computation and intentionally omits plotting, printing, interactive R interfaces, S3 presentation behavior, and other presentation-only code.

The translated package uses a single `real64` kind, exported as `dp`, and has no external numerical dependency in its current scope. It therefore does not require system BLAS/LAPACK libraries and does not vendor translated dependencies. Dense eigensystem, regression, and small-matrix operations use package-local Fortran numerical support. Existing top-level translations and shared modules in `Fortran-from-R-packages` were reviewed before implementation; no external package API was required by the retained implementation paths.

## Build and test

```text
fpm build
fpm test
fpm run --example dense_fpca
```

The package is standard free-form Fortran and is intended to build with `gfortran` on Windows and other FPM platforms. NaN-sensitive code uses `ieee_arithmetic`; do not build with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent assumptions.

## Public API highlights

`use fdapace` exposes the common `dp` kind together with integration/support conversion, basis construction, local smoothers, dynamic and concurrent correlation, bandwidth helpers, growth z-scores, simulation/sparsification, dense FPCA/FSVD, mean/covariance and cross-covariance estimators, derivative/variance-process analyses, Stringing and functional clustering, functional/additive/concurrent regression approximations, bootstrap FLM intervals, functional optimal design, FPC conditional quantiles, TVAM/VCAM approximations, and dense WFDA alignment.

The public Fortran names make narrowed interfaces explicit where useful—for example, R `FPCA` maps to `fpca_dense` and scalar-response R `FLM` maps to `flm_scalar_dense`. RNG-based routines use the Fortran intrinsic RNG; call `seed_rng` for deterministic Fortran runs, but streams do not match R.

## Translation coverage

Package status: **substantial**. **48 of 48 (100%)** exported computational R functions have a meaningful complete, substantial, or partial Fortran mapping. The 100% fraction measures **mapped computational functions**, not complete R compatibility: many mappings intentionally implement a numerically useful dense or core branch while omitting sparse/ragged PACE paths, R formula/S3/list construction, particular optimizer or GAM backends, R warnings/recycling, and other interface-specific behavior.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `BwNN` | `fdapace_statistics` | `bw_nn` | substantial |
| `ConvertSupport` | `fdapace_smoothing` | `convert_support`, `convert_support_functions`, `convert_support_covariance` | substantial |
| `CreateBasis` | `fdapace_basis` | `create_basis` | substantial |
| `cumtrapzRcpp` | `fdapace_smoothing` | `cumtrapz_rcpp` | complete |
| `Dyn_test` | `fdapace_statistics` | `dyn_test` | substantial |
| `DynCorr` | `fdapace_statistics` | `dyn_corr` | complete |
| `FAM` | `fdapace_models` | `fam` | partial |
| `FCCor` | `fdapace_statistics` | `fc_cor` | substantial |
| `FClust` | `fdapace_analysis` | `fclust` | partial |
| `FCReg` | `fdapace_models` | `fc_reg` | partial |
| `FLM` | `fdapace_models` | `flm_scalar_dense` | partial |
| `FLMCI` | `fdapace_models` | `flm_ci_scalar_dense` | partial |
| `FOptDes` | `fdapace_models` | `f_opt_des` | partial |
| `FPCA` | `fdapace_fpca` | `fpca_dense` | partial |
| `FPCAder` | `fdapace_analysis` | `fpca_der` | partial |
| `FPCquantile` | `fdapace_models` | `fpc_quantile` | partial |
| `FSVD` | `fdapace_fpca` | `fsvd_dense` | partial |
| `FVPA` | `fdapace_analysis` | `fvpa` | substantial |
| `GetCovSurface` | `fdapace_covariance` | `get_cov_surface` | partial |
| `GetCrCorYX` | `fdapace_covariance` | `get_cr_cor_yx` | complete |
| `GetCrCorYZ` | `fdapace_covariance` | `get_cr_cor_yz` | complete |
| `GetCrCovYX` | `fdapace_covariance` | `get_cr_cov_yx` | partial |
| `GetCrCovYZ` | `fdapace_covariance` | `get_cr_cov_yz` | partial |
| `GetMeanCI` | `fdapace_covariance` | `get_mean_ci` | partial |
| `GetMeanCurve` | `fdapace_covariance` | `get_mean_curve` | partial |
| `GetNormalisedSample` | `fdapace_fpca` | `get_normalised_sample` | partial |
| `GetNormalizedSample` | `fdapace_fpca` | `get_normalized_sample` | partial |
| `kCFC` | `fdapace_analysis` | `k_cfc` | partial |
| `Lwls1D` | `fdapace_smoothing` | `lwls1d` | substantial |
| `Lwls2D` | `fdapace_smoothing` | `lwls2d` | partial |
| `Lwls2DDeriv` | `fdapace_smoothing` | `lwls2d_deriv` | partial |
| `MakeBWtoZscore02y` | `fdapace_growth` | `make_bw_to_zscore_02y` | substantial |
| `MakeFPCAInputs` | `fdapace_simulation` | `make_fpca_inputs_dense` | partial |
| `MakeGPFunctionalData` | `fdapace_simulation` | `make_gp_functional_data` | substantial |
| `MakeHCtoZscore02y` | `fdapace_growth` | `make_hc_to_zscore_02y` | substantial |
| `MakeLNtoZscore02y` | `fdapace_growth` | `make_ln_to_zscore_02y` | substantial |
| `MakeSparseGP` | `fdapace_simulation` | `make_sparse_gp` | substantial |
| `MultiFAM` | `fdapace_models` | `multi_fam` | partial |
| `NormCurvToArea` | `fdapace_smoothing` | `norm_curv_to_area` | complete |
| `SBFitting` | `fdapace_models` | `sb_fitting` | partial |
| `SelectK` | `fdapace_fpca` | `select_k_fve`, `select_k_fixed` | partial |
| `Sparsify` | `fdapace_simulation` | `sparsify` | substantial |
| `Stringing` | `fdapace_analysis` | `stringing` | partial |
| `trapzRcpp` | `fdapace_smoothing` | `trapz_rcpp` | complete |
| `TVAM` | `fdapace_models` | `tvam` | partial |
| `VCAM` | `fdapace_models` | `vcam` | partial |
| `WFDA` | `fdapace_analysis` | `wfda` | partial |
| `Wiener` | `fdapace_simulation` | `wiener` | substantial |

There are no exported computational R functions on the untranslated list. Mapping limitations are recorded per function in `fpm.toml`, and `docs/API_COVERAGE.md` summarizes the main compatibility gaps.

## Numerical compatibility notes

The local smoothers reproduce the upstream Epanechnikov, rectangular, Gaussian, Gaussian-variance, and quartic formulas. The 2-D Gaussian-variance kernel intentionally preserves fdapace's different first- and second-coordinate variance factors. `MakeLNtoZscore02y` preserves the literal upstream male reference value `6.0477` at month 22 rather than silently correcting it. The three growth helpers accept the documented age range 0 through 24 months, avoiding the R source's `time[24]` one-based-index boundary artifact.

Dense `FPCA` covers common-grid cross-sectional covariance, optional second-difference error-variance estimation, nonnegative eigendecomposition, trapezoidal eigenfunction normalization, FVE selection, integration scores, optional shrinkage, fitted covariance/correlation, and fitted curves. Dense `FSVD` covers paired regular samples. Higher-level mappings such as `FAM`, `FLM`, `FClust`, `kCFC`, `TVAM`, `VCAM`, and `WFDA` provide deterministic dense/core numerical paths rather than reproducing every R backend or sparse-data option.

## Provenance and licensing

The upstream metadata and R/C++ source files relevant to mapped functions are retained under `upstream/` for provenance. See `NOTICE` and `LICENSE`. No dependency source is copied into this package.
