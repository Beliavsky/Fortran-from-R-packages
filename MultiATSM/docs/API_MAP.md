# API coverage map

This map groups the R entry points by computational role. Several R functions
mainly unpack lists, attach labels, perform input checks, or dispatch to the
same numerical kernel; these map to one matrix-first Fortran routine.

## Direct native counterparts

| R functionality | Fortran counterpart |
|---|---|
| `pca_weights_one_country` | `pca_weights_one_country` |
| `Spanned_Factors` | `spanned_factors` |
| `VarianceExplained` | `pca_variance_explained` |
| `VAR`, `Est_RestOLS` | `fit_var`, `fit_restricted_ols` |
| `VARX` | `fit_varx_system` |
| `BuildGVAR`, `GVAR` numerical core | `build_gvar`, `fit_gvar` |
| `Transition_Matrix` | `transition_matrix_year`, `transition_matrix_mean` |
| foreign/star-factor construction | `build_star_factors` |
| JLL orthogonalization and dynamics | `fit_jll` |
| `FeedbackMatrixRestrictionsJLL` | `jll_feedback_restrictions` |
| `CholRestrictionsJLL` | `jll_cholesky_mask` |
| `Compute_BnX_AnX`, `Get__BnXAnX` | `affine_yield_loadings`, `multicountry_affine_loadings` |
| `Rotate_Lat_Obs` | `rotate_latent_to_observed`, `pricing_factor_loadings` |
| `Get_r0` numerical core | `estimate_long_run_short_rate` |
| `MLEdensity`, `GaussianDensity`, `Get_llk` | `atsm_log_likelihood`, `gaussian_log_density` |
| `Get_SigmaYields` | `yield_error_variance` |
| `YieldsFit`, `Y_Fit` | `fitted_yields` |
| `ForecastYields`, `Gen_Forecast_Yields`, `RMSE` | `forecast_yields`, `forecast_rmse` |
| `ComputeIRFs`, `ComputeGIRFs` | `impulse_responses`, `generalized_impulse_responses` |
| `ComputeFEVDs`, `ComputeGFEVDs` | `forecast_error_variance_decomposition`, `generalized_fevd` |
| `ExpectedComponent`, `Compute_EP` | `expected_short_rate_component` |
| `TermPremia`, `ForwardPremia` | `term_premium`, `forward_rates` |
| `OptimizationSetup_ATSM` | `bfgs_minimize`, `nelder_mead_minimize` |
| `Jac_approx`, `richardson_diff` | `numerical_jacobian`, `numerical_gradient` |
| `scaling_from_jacobian` | `scale_from_jacobian` |
| stationary Jordan constraint purpose | `stabilize_transition` |
| `Aux_PSD`, `True_PSD` | `psd_to_lower_parameters`, `lower_factor_to_psd` |
| `Aux_BlockDiag`, `True_BlockDiag` | `block_diagonal_psd` |
| `ResampleResiduals_BS` | `resample_residuals` |
| artificial VAR-series generation | `simulate_var` |
| VAR bootstrap parameter draws | `bootstrap_var` |
| bootstrap quantile bounds | `percentile_bounds` |
| `Bias_Correc_VAR`, `SA_algorithm`, `Gen_art_series` | `bias_correct_var` |
| `shrink_FeedMat_BC` | `shrink_transition` |

## Adapted interfaces

### Model assembly

The R package represents each model with deeply nested lists and character
labels. The Fortran port uses derived types (`var_model`, `gvar_model`,
`jll_model`, `affine_loadings`) and explicit matrices. The numerical equations
are available, but there is no string-based `ModelType` dispatcher equivalent
to `Optimization`, `NumOutputs`, or `Bootstrap`.

### Parameter constraints

The R code maps full matrices to auxiliary Jordan-eigenvalue vectors and back.
The Fortran port exposes the numerically important constraints directly:
transition matrices can be projected to a chosen spectral-radius bound, and
covariances are built from lower-triangular factors. This is more transparent
for Fortran callers but is not an identical auxiliary-vector ordering.

### JLL covariance restrictions

`fit_jll` reproduces the factor orthogonalization, feedback restrictions,
PI transformation, and covariance masking. The upstream R code can perform an
inner likelihood optimization of the restricted Cholesky entries. The Fortran
routine currently projects the unrestricted Cholesky factor onto the JLL mask;
it does not repeat that inner Nelder-Mead optimization.

### Multi-country affine loadings

`multicountry_affine_loadings` retains full cross-country B loadings. This is a
more general numerical representation than helper paths in the R code that
extract only country-specific blocks for selected outputs.

### Bias correction

`bias_correct_var` implements the Bauer-Rudebusch-Wu-style stochastic
approximation for an unrestricted matrix VAR(1). R dispatch that re-estimates a
GVAR or JLL restriction structure inside every draw remains the caller's
responsibility: simulate a draw, call `fit_gvar` or `fit_jll`, and feed the
result into a custom iteration.

## Deliberately omitted

The following are R/application infrastructure rather than reusable Fortran
numerical kernels:

- S3 classes, `print`, `summary`, `plot`, and `autoplot` methods;
- ggplot/cowplot graph construction and filesystem folder creation;
- R list, character-label, row-name, and formula dispatch;
- Excel import, `.rda` datasets, CRAN examples, date parsing, and frequency
  labels;
- data-frame merging, missing-value cleaning, and interactive input checks;
- package-specific HTML/bookdown/vignette generation;
- confidence-bound graph formatting (numeric quantiles are provided);
- model-output serialization tied to R object layouts.
