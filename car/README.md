# car

Modern free-form Fortran translation of computational kernels from the R package **car** (Companion to Applied Regression), packaged for the Fortran Package Manager (FPM).

This translation is based on upstream **car 3.1-5**. It focuses on reusable numerical algorithms and deliberately omits plotting, interactive graphics, formula parsing, R object construction, and most S3 model dispatch. The original R sources used for provenance are retained under `upstream/`.

## Build and test

From the `car` directory:

```text
fpm build
fpm test
fpm run --example basic
```

The package has no external dependencies and does not link system BLAS or LAPACK. It therefore remains portable to Windows gfortran without separately installed numerical libraries.

## Implemented areas

The public `car` module re-exports the translated API. Major implemented areas include Box-Cox, Yeo-Johnson, and BCN power transforms; contrast matrices; adaptive kernel density estimation; ellipse coordinate generation; weighted cross-products; HC0-HC4 covariance estimators; GVIFs; numeric linear-hypothesis tests; delta-method standard errors; Durbin-Watson statistics; Levene and NCV tests; Bonferroni outlier probabilities; numeric recoding and grouped reductions; Box-Cox profiling; power-parameter estimation; and inverse-transformation estimation.

The translation uses one real kind, `dp`, defined once from `iso_fortran_env::real64` in `car_kinds` and re-exported from `car`. No maintained Fortran source uses `double precision`, `real*8`, D-exponent literals, or implicit interfaces.

## Quick example

```fortran
program example
   use car, only: bc_power, dp
   implicit none
   real(dp) :: y(3), z(3)
   integer :: status

   y = [1.0_dp, 2.0_dp, 4.0_dp]
   call bc_power(y, 0.0_dp, z, status)
   if (status /= 0) error stop "bc_power failed"
   print *, z
end program example
```

Most routines report invalid inputs through integer `status` outputs rather than R conditions. NaN-aware translated behavior uses `ieee_is_nan` from `ieee_arithmetic`.

## Translation coverage

Package status: **substantial**.

**51 of 51 (100.0%)** distinct exported computational R functions have a meaningful complete, substantial, or partial Fortran mapping. This fraction measures mapped computational functions, **not** complete compatibility with R's formulas, S3 classes, object semantics, warnings, graphics, RNG, or model ecosystem.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `contr.Helmert` | `car` | `helmert_contrast` | substantial |
| `contr.Sum` | `car` | `sum_contrast` | substantial |
| `contr.Treatment` | `car` | `treatment_contrast` | substantial |
| `ellipse` | `car` | `ellipse_coordinates` | substantial |
| `logit` | `car` | `logit_transform` | substantial |
| `recode` | `car` | `recode_numeric` | partial |
| `Recode` | `car` | `recode_numeric` | partial |
| `sigmaHat` | `car` | `sigma_hat` | partial |
| `basicPower` | `car` | `basic_power` | substantial |
| `bcPower` | `car` | `bc_power` | substantial |
| `yjPower` | `car` | `yj_power` | substantial |
| `powerTransform` | `car` | `power_transform_fit` | partial |
| `bcnPower` | `car` | `bcn_power` | substantial |
| `bcnPowerInverse` | `car` | `bcn_power_inverse` | complete |
| `boxCox` | `car` | `box_cox_profile` | substantial |
| `boxCoxVariable` | `car` | `box_cox_variable` | complete |
| `dbiwt` | `car` | `biweight_kernel` | complete |
| `deltaMethod` | `car` | `delta_method_gradient` | partial |
| `depan` | `car` | `epanechnikov_kernel` | complete |
| `dwt` | `car` | `durbin_watson` | partial |
| `durbinWatsonTest` | `car` | `durbin_watson` | partial |
| `testTransform` | `car` | `power_transform_test` | partial |
| `hccm` | `car` | `hccm_matrix` | substantial |
| `invTranEstimate` | `car` | `inv_transform_estimate` | substantial |
| `leveneTest` | `car` | `levene_test` | substantial |
| `lht` | `car` | `linear_hypothesis_test` | partial |
| `linearHypothesis` | `car` | `linear_hypothesis_test` | partial |
| `ncvTest` | `car` | `ncv_test` | partial |
| `outlierTest` | `car` | `outlier_bonferroni` | partial |
| `vif` | `car` | `gvif_from_covariance` | substantial |
| `adaptiveKernel` | `car` | `adaptive_kernel` | substantial |
| `Tapply` | `car` | `tapply_numeric` | partial |
| `strings2factors` | `car` | `strings2factors_codes` | partial |
| `some` | `car` | `some_indices` | partial |
| `which.names` | `car` | `which_names` | partial |
| `whichNames` | `car` | `which_names` | partial |
| `compareCoefs` | `car` | `compare_coefs_numeric` | partial |
| `linearHypothesis.default` | `car` | `linear_hypothesis_test` | substantial |
| `wcrossprod` | `car` | `wcrossprod` | substantial |
| `boxTidwell` | `car` | `box_tidwell` | substantial |
| `poTest` | `car` | `po_test` | substantial |
| `Predict` | `car` | `predict_lm` | substantial |
| `Confint` | `car` | `confint_lm` | partial |
| `S` | `car` | `lm_summary` | partial |
| `Anova` | `car` | `anova_lm_type3` | partial |
| `makeHypothesis` | `car` | `make_hypothesis` | partial |
| `Manova` | `car` | `manova_type3_ssp`, `manova_test_statistics` | partial |
| `bootCase` | `car` | `lm_bootstrap_case` | partial |
| `nextBoot` | `car` | `lm_refit_rows` | partial |
| `matchCoefs` | `car` | `match_coefficient_names` | partial |
| `Boot` | `car` | `lm_bootstrap_case` | partial |

All functions in the stated coverage basis now have a mapping, but many mappings intentionally cover only reusable numerical kernels. `fpm.toml` is the machine-readable source of truth for the function-by-function coverage map. `docs/API_COVERAGE.md` explains the coverage basis and compatibility limits in more detail.

## Important compatibility limits

The Fortran API generally accepts numerical arrays that correspond to components R would extract from fitted model objects. For example, `hccm_matrix` accepts a model matrix, residuals, unscaled covariance matrix, and leverage values; `linear_hypothesis_test` accepts coefficients, their covariance matrix, and a numeric hypothesis matrix. This avoids reproducing R formulas, environments, S3 dispatch, and model classes.

`power_transform_fit` currently estimates a single response transformation parameter, while R `powerTransform` also supports multivariate responses and richer Hessian-based inference. `box_tidwell` translates the numeric default iteration, initial score tests, and omnibus test but not formula/model-frame handling. `po_test` accepts an explicit design matrix and ordered response rather than a fitted `polr` object. `durbin_watson` implements deterministic residual-vector statistics but not the simulation-based p-values of the `lm` method. `recode_numeric` handles numeric interval rules rather than R expression parsing and factor construction.

The data-oriented mappings are also intentionally narrower than R: `strings2factors_codes` operates on homogeneous character matrices rather than heterogeneous data frames, `some_indices` returns sampled positions and does not reproduce R RNG streams, `which_names` performs exact lookup on explicit names, and `compare_coefs_numeric` accepts aligned coefficient/covariance arrays rather than fitted model objects.

The newer model-level mappings remain intentionally narrower than their R generics: `lm_summary`, `confint_lm`, and `predict_lm` operate on explicit numeric linear-model components; `anova_lm_type3` implements Type III LM tests; `manova_type3_ssp` plus `manova_test_statistics` implement Type III SSP construction and Pillai/Wilks/Hotelling-Lawley tests but not Type II, repeated measures, or Roy's root; `lm_bootstrap_case` requires caller-supplied bootstrap indices rather than owning RNG/parallel model updates; and `make_hypothesis` accepts ordinary linear-expression syntax but not every nonstandard R coefficient name. See the coverage notes in `fpm.toml` for function-specific limitations.

## Repository integration

The target repository is `Beliavsky/Fortran-from-R-packages`. Existing shared facilities and top-level translations were checked before implementation. The translated slice here is intentionally dependency-free: none of its implemented public algorithms requires the higher-level `rfortran-linalg` API, direct LAPACK, ARPACK, or an existing translated R-package dependency. No BLAS, LAPACK, ARPACK, `r.f90`, `r_mod.f90`, or other package sources are vendored.

## Provenance and license

Upstream `car` is authored by John Fox, Sanford Weisberg, Brad Price, and numerous contributors and is distributed under **GPL (>= 2)**. The Fortran translation is distributed under the same GPL-2.0-or-later terms. See `NOTICE.md`, `LICENSE`, `upstream/DESCRIPTION`, and `upstream/inst/CITATION`.
