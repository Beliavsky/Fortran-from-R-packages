# API coverage

## Coverage basis

The translation uses the required basis **exported computational R functions**. The denominator is 51 distinct exported functions from upstream `NAMESPACE` that perform statistical, numerical, fitting, sampling, name-selection, or substantive data-transformation operations. Defunct aliases are excluded because they only raise `.Defunct`; plotting/graphics, printing/formatting, browser/data-download, interactive-display, and presentation-only exports are excluded. Utilities such as `strings2factors`, `some`, `whichNames`, `which.names`, and `compareCoefs` are included because they return transformed, selected, or statistically derived data rather than merely formatting output.

The package status is **substantial**. All 51 functions in the coverage basis now have a meaningful mapping, but several central model-level facilities remain partial because the Fortran API works with explicit numerical arrays rather than reproducing R formula/S3/model-object semantics.

## Summary

- R functions in coverage basis: **51**
- R functions mapped: **51**
- Fraction mapped: **1.0 (100.0%)**
- Untranslated: none in the stated coverage basis

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

## Public Fortran modules

- `car`: facade module and recommended import surface; re-exports `dp` and all translated public procedures.
- `car_transforms`: power transformations, inverse BCN transformation, logit, and Box-Cox variable.
- `car_power_fit`: univariate power fitting, transformation tests, Box-Cox profiles, and inverse-transformation estimation.
- `car_box_tidwell`: numeric Box-Tidwell transformation fitting, score tests, and omnibus testing.
- `car_po_test`: proportional-odds diagnostics using threshold logistic fits.
- `car_lm`: explicit-array linear-model summaries, coefficient intervals, predictions, case bootstrap, and sampled-row refits.
- `car_anova`: Type III linear-model term tests.
- `car_manova`: Type III multivariate hypothesis/error SSP matrices and Pillai/Wilks/Hotelling-Lawley test statistics.
- `car_hypothesis`: constrained text-to-numeric linear-hypothesis parsing.
- `car_contrasts`: treatment, sum, and Helmert contrast matrices.
- `car_kernels`: normal, Epanechnikov, biweight, normal-reference bandwidth, and adaptive density estimation.
- `car_geometry`: ellipse coordinate generation.
- `car_matrix`: weighted cross-products, residual scale, HC covariance, GVIF, numeric Wald tests, and delta-method variance.
- `car_diagnostics`: Durbin-Watson, Levene, NCV, and outlier probability kernels.
- `car_data`: numeric recoding, grouped reductions, name lookup, sampling indices, character-to-factor coding, and coefficient-comparison statistics.
- `car_utils`: package-private numerical support for probability tails, quantiles, simple linear algebra, and optimization.

## Interpretation of mapping status

`complete` means the numerical operation itself is translated without a material computational omission for the stated R function. `substantial` means the principal numerical operation is present but R object/model extraction or secondary options are omitted. `partial` means a useful numerical core is present while important supported R modes remain absent. A complete function mapping does not imply R interface compatibility.
