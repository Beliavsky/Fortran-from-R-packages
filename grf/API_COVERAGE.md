# API coverage

## Coverage basis

The coverage denominator follows the project convention: distinct exported computational R functions, plus explicitly exported computational S3 methods. Plotting/printing/presentation routines, the dispatch-only `get_scores` generic, and configuration-only helpers are excluded.

**Package status: substantial — 41 of 41 (100.0%) mapped computational functions.**

The percentage measures whether a meaningful Fortran implementation is mapped to each counted R operation; it does **not** claim complete R compatibility. Current major translated parity features include OOB nuisance forests, grouped little-bag sampling and debiased variances for the supported scalar families, causal-survival censoring/IPCW pseudo-outcomes, IV compliance forests for binary-instrument scores, boosted OOB stopping, regression hyperparameter tuning, local-linear ridge tuning/residual splitting, cluster-aware sampling and inference, all twelve causal simulation selectors, and all six causal-survival simulation selectors.

The principal remaining computational gaps are family-specific tuning outside regression/local-linear/boosted regression; exact alternative sandwich/TMLE/subset pathways; some specialized survival and causal-survival inference behavior; exact RATE weighting/bootstrap variants; exact multi-arm covariance variants; and exact R RNG/threading behavior. The per-function notes in `fpm.toml` are authoritative for mapped differences.

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
