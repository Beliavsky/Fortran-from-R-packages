# grf — modern Fortran translation

This directory translates the computational core of the R package **grf 2.6.1** (*Generalized Random Forests*) to modern free-form Fortran with FPM. The implementation is self-contained: it does not require R, Rcpp, RcppEigen, BLAS, LAPACK, or a separately installed numerical library.

The common tree engine supports randomized `mtry`, subsampling, honesty, missing-value routing, balance constraints, grouped "little-bag" tree sampling, family-specific node relabeling, honest leaf repopulation/pruning, adaptive forest weights, and local moment prediction. Public APIs cover regression, multi-response regression, probability, quantile, survival, causal, instrumental, causal-survival, multi-treatment/local-coefficient, local-linear regression, and boosted regression workflows, together with analysis and simulation helpers.

## Build

From this top-level directory:

```text
fpm build
fpm test
```

Examples are discovered by FPM from `example/`. All maintained Fortran source is free-form and uses the single public real kind `dp`.

## Major parity features

The current translation includes several mechanisms that were partial in the first release:

- automatic OOB nuisance forests for causal, IV, multi-arm, and local-coefficient estimators;
- GRF-style grouped half-sample / little-bag training through `ci_group_size`, including tree-count rounding and objective-Bayes variance debiasing for regression, causal, IV, causal-survival, and local-linear predictions;
- censoring/event/treatment nuisance survival forests and IPCW/doubly robust pseudo-outcomes for causal-survival RMST and survival-probability targets;
- conditional treatment-variance causal scores and binary-instrument compliance-forest IV scores;
- boosted-regression OOB debiased-error stopping and reuse of first-stage tuned forest controls;
- deterministic regression hyperparameter tuning over the upstream seven-parameter ranges with a pure-Fortran Gaussian-process surrogate;
- local-linear OOB ridge-path selection, covariance-weighted ridge penalties, and optional training-time local-linear residual splitting with split lambda/variables/cutoff controls; and
- all twelve upstream causal simulation benchmark designs, including upstream-style baseline, treatment-effect, and noise scaling; and
- all six upstream causal-survival simulation designs, including supplied/correlated predictors and deterministic Monte Carlo target effects.

## Important compatibility notes

This is a numerical translation, not an R compatibility layer. R lists/S3 dispatch, factor labels, formula/data-frame processing, printing/plotting, C++ threading, and exact R RNG streams are intentionally omitted. Family-specific automatic hyperparameter tuning is currently implemented for regression (and reused by local-linear/boosted regression where applicable), not for every GRF family. Several alternative sandwich/TMLE pathways, exact RATE bootstrap details, and some specialized survival prediction/inference options remain partial. The pure-Fortran tuning surrogate follows GRF's parameter ranges and OOB error selection but is not an exact reimplementation of R's `DiceKriging` optimization.

## Translation coverage

**Status: substantial — 41 of 41 (100.0%) mapped computational functions.**

The 100.0% fraction means every function in the stated coverage basis has a meaningful complete/substantial/partial Fortran mapping. It does **not** mean complete R interface or estimator-option compatibility. Of the 41 mappings, 33 are currently marked substantial and 8 partial. See `API_COVERAGE.md` and the per-function `[[extra.translation.function]]` records in `fpm.toml` for material differences.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `average_treatment_effect` | `grf_analysis` | `average_treatment_effect` | partial |
| `best_linear_projection` | `grf_analysis` | `best_linear_projection` | partial |
| `boosted_regression_forest` | `grf_train` | `boosted_regression_forest` | substantial |
| `causal_forest` | `grf_train` | `causal_forest` | substantial |
| `causal_survival_forest` | `grf_train` | `causal_survival_forest` | substantial |
| `generate_causal_data` | `grf_simulation` | `generate_causal_data` | substantial |
| `generate_causal_survival_data` | `grf_simulation` | `generate_causal_survival_data` | substantial |
| `get_forest_weights` | `grf_analysis` | `get_forest_weights` | substantial |
| `get_leaf_node` | `grf_analysis` | `get_leaf_node` | substantial |
| `get_tree` | `grf_analysis` | `get_tree` | substantial |
| `instrumental_forest` | `grf_train` | `instrumental_forest` | substantial |
| `ll_regression_forest` | `grf_train` | `ll_regression_forest` | substantial |
| `lm_forest` | `grf_train` | `lm_forest` | substantial |
| `merge_forests` | `grf_analysis` | `merge_forests` | substantial |
| `multi_arm_causal_forest` | `grf_train` | `multi_arm_causal_forest` | substantial |
| `multi_regression_forest` | `grf_train` | `multi_regression_forest` | substantial |
| `probability_forest` | `grf_train` | `probability_forest` | substantial |
| `quantile_forest` | `grf_train` | `quantile_forest` | substantial |
| `rank_average_treatment_effect` | `grf_analysis` | `rank_average_treatment_effect` | partial |
| `rank_average_treatment_effect.fit` | `grf_analysis` | `rank_average_treatment_effect_fit` | partial |
| `regression_forest` | `grf_train` | `regression_forest` | substantial |
| `split_frequencies` | `grf_analysis` | `split_frequencies` | substantial |
| `survival_forest` | `grf_train` | `survival_forest` | partial |
| `test_calibration` | `grf_analysis` | `test_calibration` | partial |
| `variable_importance` | `grf_analysis` | `variable_importance` | substantial |
| `get_scores.causal_forest` | `grf_analysis` | `get_scores_causal_forest` | substantial |
| `get_scores.causal_survival_forest` | `grf_analysis` | `get_scores_causal_survival_forest` | substantial |
| `get_scores.instrumental_forest` | `grf_analysis` | `get_scores_instrumental_forest` | substantial |
| `get_scores.multi_arm_causal_forest` | `grf_analysis` | `get_scores_multi_arm_causal_forest` | partial |
| `predict.boosted_regression_forest` | `grf_predict` | `predict_boosted_regression_forest` | substantial |
| `predict.causal_forest` | `grf_predict` | `predict_causal_forest` | substantial |
| `predict.causal_survival_forest` | `grf_predict` | `predict_causal_survival_forest` | substantial |
| `predict.instrumental_forest` | `grf_predict` | `predict_instrumental_forest` | substantial |
| `predict.ll_regression_forest` | `grf_predict` | `predict_ll_regression_forest` | substantial |
| `predict.lm_forest` | `grf_predict` | `predict_lm_forest` | substantial |
| `predict.multi_arm_causal_forest` | `grf_predict` | `predict_multi_arm_causal_forest` | substantial |
| `predict.multi_regression_forest` | `grf_predict` | `predict_multi_regression_forest` | substantial |
| `predict.probability_forest` | `grf_predict` | `predict_probability_forest` | substantial |
| `predict.quantile_forest` | `grf_predict` | `predict_quantile_forest` | substantial |
| `predict.regression_forest` | `grf_predict` | `predict_regression_forest` | substantial |
| `predict.survival_forest` | `grf_predict` | `predict_survival_forest` | partial |

## Provenance and license

The supplied upstream package is retained under `upstream/`, except for the generated upstream build artifact documented in `VALIDATION.md`. See `UPSTREAM.md` and `NOTICE.md`. The package metadata declares GPL-3; this translation is distributed under GPL-3.0-only consistently with that metadata.
