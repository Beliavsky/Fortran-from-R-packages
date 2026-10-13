# Computational API coverage audit

Source: CRAN loo 2.10.1 NAMESPACE. Package status: substantial.
There are 65 explicit export() names, 21 documented exclusions, 44 computational
functions, 40 mapped (90.91%), and four untranslated. Names registered only with
S3method() do not add exported-function denominator entries. Every exported
computational method is counted separately when explicitly exported.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `.compute_point_estimate` | `loo_evaluation` | `compute_point_estimate` | complete |
| `.thin_draws` | `loo_evaluation` | `thin_draws` | substantial |
| `E_loo` | `loo_diagnostics` | `E_loo` | substantial |
| `compare` | `loo_diagnostics` | `compare` | partial |
| `crps` | `loo_metrics` | `crps` | substantial |
| `elpd` | `loo_evaluation` | `elpd` | substantial |
| `kfold_split_grouped` | `loo_folds` | `kfold_split_grouped` | substantial |
| `kfold_split_random` | `loo_folds` | `kfold_split_random` | substantial |
| `kfold_split_stratified` | `loo_folds` | `kfold_split_stratified` | substantial |
| `loo` | `loo_evaluation` | `loo` | substantial |
| `loo.array` | `loo_evaluation` | `loo_array` | substantial |
| `loo.function` | `loo_evaluation` | `loo_function` | substantial |
| `loo.matrix` | `loo_evaluation` | `loo_matrix` | substantial |
| `loo_approximate_posterior` | `loo_evaluation` | `loo_approximate_posterior` | substantial |
| `loo_approximate_posterior.array` | `loo_evaluation` | `loo_approximate_posterior_array` | substantial |
| `loo_approximate_posterior.function` | `loo_evaluation` | `loo_approximate_posterior_function` | substantial |
| `loo_approximate_posterior.matrix` | `loo_evaluation` | `loo_approximate_posterior_matrix` | substantial |
| `loo_compare` | `loo_diagnostics` | `loo_compare` | substantial |
| `loo_crps` | `loo_metrics` | `loo_crps` | substantial |
| `loo_i` | `loo_evaluation` | `loo_i` | substantial |
| `loo_model_weights` | `loo_model_weights` | `model_weights` | partial |
| `loo_model_weights.default` | `loo_model_weights` | `model_weights` | partial |
| `loo_predictive_metric` | `loo_metrics` | `loo_predictive_metric` | substantial |
| `loo_scrps` | `loo_metrics` | `loo_scrps` | substantial |
| `mcse_loo` | `loo_diagnostics` | `mcse_loo` | complete |
| `pareto_k_ids` | `loo_diagnostics` | `pareto_k_ids` | complete |
| `pareto_k_table` | `loo_diagnostics` | `pareto_k_table` | substantial |
| `pseudobma_weights` | `loo_model_weights` | `pseudobma_weights` | substantial |
| `psis` | `loo_importance` | `psis` | substantial |
| `psislw` | `loo_importance` | `psislw` | substantial |
| `relative_eff` | `loo_evaluation` | `relative_eff`, `relative_eff_function` | substantial |
| `scrps` | `loo_metrics` | `scrps` | substantial |
| `sis` | `loo_importance` | `sis` | substantial |
| `stacking_weights` | `loo_model_weights` | `stacking_weights` | substantial |
| `tis` | `loo_importance` | `tis` | substantial |
| `waic` | `loo_evaluation` | `waic` | substantial |
| `waic.array` | `loo_evaluation` | `waic_array` | substantial |
| `waic.function` | `loo_evaluation` | `waic_function` | substantial |
| `waic.matrix` | `loo_evaluation` | `waic_matrix` | substantial |
| `weights.importance_sampling` | `loo_importance` | `importance_weights` | complete |

## Mapping notes

- `.compute_point_estimate`: Column means of a numeric posterior-draw matrix; the single R output row is represented by a vector.
- `.thin_draws`: Deterministic draw subsampling preserves R integer-division precedence; arbitrary S3 objects and NULL pass-through are omitted.
- `E_loo`: Weighted mean, corrected variance, SD, type-seven/uniform or interpolated weighted quantiles, and function-specific Pareto k; explicit PSIS result and optional raw ratios replace S3 arguments. Centered variance improves release numerical stability.
- `compare`: Deprecated two-model pointwise difference and paired standard error are implemented. The old multi-model table interface is omitted; use loo_compare for ordered multi-model comparisons.
- `crps`: Numeric predictive score and standard error. Explicit permutation indexes replace R RNG; LOO variants take normalized marginal and separately PSIS-smoothed paired joint weights. R object construction and automatic weight fitting are omitted.
- `elpd`: Matrix and iteration/chain-array log predictive density, information criterion, totals and standard errors; inputs are finite and R S3 metadata are omitted.
- `kfold_split_grouped`: Fold assignment algorithm with explicitly supplied permutations in place of R random-number generation. Categories/groups are integer labels, ordered numerically; invalid arguments return zeros.
- `kfold_split_random`: Fold assignment algorithm with explicitly supplied permutations in place of R random-number generation. Categories/groups are integer labels, ordered numerically; invalid arguments return zeros.
- `kfold_split_stratified`: Fold assignment algorithm with explicitly supplied permutations in place of R random-number generation. Categories/groups are integer labels, ordered numerically; invalid arguments return zeros.
- `loo`: PSIS/SIS/TIS pointwise LOO, effective parameter count, information criterion, MCSE and summary errors. Arrays are iterations by chains by observations; callbacks are explicit procedure pointers with indexed draw-vector output. All likelihoods must be finite; R objects, parallel dispatch and warnings are omitted.
- `loo.array`: PSIS/SIS/TIS pointwise LOO, effective parameter count, information criterion, MCSE and summary errors. Arrays are iterations by chains by observations; callbacks are explicit procedure pointers with indexed draw-vector output. All likelihoods must be finite; R objects, parallel dispatch and warnings are omitted.
- `loo.function`: PSIS/SIS/TIS pointwise LOO, effective parameter count, information criterion, MCSE and summary errors. Arrays are iterations by chains by observations; callbacks are explicit procedure pointers with indexed draw-vector output. All likelihoods must be finite; R objects, parallel dispatch and warnings are omitted.
- `loo.matrix`: PSIS/SIS/TIS pointwise LOO, effective parameter count, information criterion, MCSE and summary errors. Arrays are iterations by chains by observations; callbacks are explicit procedure pointers with indexed draw-vector output. All likelihoods must be finite; R objects, parallel dispatch and warnings are omitted.
- `loo_approximate_posterior`: PSIS ratios corrected by target-minus-approximation log densities, with full LOO summaries. Density vectors must be finite and use the same draw order; R dispatch, deprecated log_q and class metadata are omitted.
- `loo_approximate_posterior.array`: PSIS ratios corrected by target-minus-approximation log densities, with full LOO summaries. Density vectors must be finite and use the same draw order; R dispatch, deprecated log_q and class metadata are omitted.
- `loo_approximate_posterior.function`: PSIS ratios corrected by target-minus-approximation log densities, with full LOO summaries. Density vectors must be finite and use the same draw order; R dispatch, deprecated log_q and class metadata are omitted.
- `loo_approximate_posterior.matrix`: PSIS ratios corrected by target-minus-approximation log densities, with full LOO summaries. Density vectors must be finite and use the same draw order; R dispatch, deprecated log_q and class metadata are omitted.
- `loo_compare`: Orders at least two models by total predictive score; returns paired differences/SEs, normal-approximation worse probabilities and each model ELPD/SE. Takes pointwise scores rather than R objects; ancillary penalty/IC columns, diagnostic text, subsampling comparison and many-model selection-bias warning are omitted.
- `loo_crps`: Numeric predictive score and standard error. Explicit permutation indexes replace R RNG; LOO variants take normalized marginal and separately PSIS-smoothed paired joint weights. R object construction and automatic weight fitting are omitted.
- `loo_i`: PSIS/SIS/TIS pointwise LOO, effective parameter count, information criterion, MCSE and summary errors. Arrays are iterations by chains by observations; callbacks are explicit procedure pointers with indexed draw-vector output. All likelihoods must be finite; R objects, parallel dispatch and warnings are omitted.
- `loo_model_weights`: Routes stacking, pseudo-BMA and pseudo-BMA+ for supplied pointwise log predictive scores. The R list-of-model/log-likelihood evaluation layer is omitted; Bayesian-bootstrap weights must be supplied explicitly.
- `loo_model_weights.default`: Routes stacking, pseudo-BMA and pseudo-BMA+ for supplied pointwise log predictive scores. The R list-of-model/log-likelihood evaluation layer is omitted; Bayesian-bootstrap weights must be supplied explicitly.
- `loo_predictive_metric`: MAE, MSE, RMSE, accuracy and balanced accuracy with mean-based PSIS-weighted predictions. Weights are supplied explicitly; optional upstream E_loo quantile/variance prediction selection is omitted.
- `loo_scrps`: Numeric predictive score and standard error. Explicit permutation indexes replace R RNG; LOO variants take normalized marginal and separately PSIS-smoothed paired joint weights. R object construction and automatic weight fitting are omitted.
- `mcse_loo`: Quadrature aggregation of pointwise LOO Monte Carlo errors gated by the sample-size-dependent or supplied Pareto-k threshold; unreliable estimates return quiet NaN.
- `pareto_k_ids`: One-based indexes above the supplied or automatic reliability threshold, operating on explicit importance results.
- `pareto_k_table`: Three numeric reliability groups with counts, proportions and minimum reliable ESS. Human-readable row/column names and formatted printing are omitted.
- `pseudobma_weights`: Plain pseudo-BMA and Bayesian-bootstrap pseudo-BMA+ weighting; explicit normalized bootstrap rows replace the R Dirichlet RNG and B/alpha controls.
- `psis`: Vector, matrix and chain-array importance sampling with normalized/raw log weights, relative-efficiency-corrected ESS and tail diagnostics. PSIS uses release ceiling tail length and sibling posterior GPD fitting; warnings and parallel R dispatch are omitted.
- `psislw`: Deprecated quantile-tail PSIS algorithm including weak-prior-off GPD fit, selectable tail proportion and truncation exponent; normalized output, matrix input. Deprecated llfun/llargs/internal COMPUTE_LOOS modes are omitted.
- `relative_eff`: Likelihood-scale ESS divided by draw count via sibling posterior ess_mean. Vector/matrix chain IDs preserve within-chain order even when interleaved; arrays and explicit callbacks supported. Equal-length contiguous positive chain IDs required; parallel R dispatch and object getter omitted.
- `scrps`: Numeric predictive score and standard error. Explicit permutation indexes replace R RNG; LOO variants take normalized marginal and separately PSIS-smoothed paired joint weights. R object construction and automatic weight fitting are omitted.
- `sis`: Vector, matrix and chain-array importance sampling with normalized/raw log weights, relative-efficiency-corrected ESS and tail diagnostics. PSIS uses release ceiling tail length and sibling posterior GPD fitting; warnings and parallel R dispatch are omitted.
- `stacking_weights`: Same concave log-score objective on the closed probability simplex, solved by vertex exchange and exact line search instead of constrOptim. Reports a primal-dual gap and explicit nonconvergence; R optimizer-specific controls are omitted.
- `tis`: Vector, matrix and chain-array importance sampling with normalized/raw log weights, relative-efficiency-corrected ESS and tail diagnostics. PSIS uses release ceiling tail length and sibling posterior GPD fitting; warnings and parallel R dispatch are omitted.
- `waic`: WAIC pointwise ELPD, variance penalty, information criterion and summary errors for finite log likelihoods with at least two draws. Matrix/chain-array or indexed callback APIs replace R object/formal-argument handling; warnings omitted.
- `waic.array`: WAIC pointwise ELPD, variance penalty, information criterion and summary errors for finite log likelihoods with at least two draws. Matrix/chain-array or indexed callback APIs replace R object/formal-argument handling; warnings omitted.
- `waic.function`: WAIC pointwise ELPD, variance penalty, information criterion and summary errors for finite log likelihoods with at least two draws. Matrix/chain-array or indexed callback APIs replace R object/formal-argument handling; warnings omitted.
- `waic.matrix`: WAIC pointwise ELPD, variance penalty, information criterion and summary errors for finite log likelihoods with at least two draws. Matrix/chain-array or indexed callback APIs replace R object/formal-argument handling; warnings omitted.
- `weights.importance_sampling`: Selects normalized or raw weights on the log or ordinary scale. Centered stored normalized log weights preserve mass under enormous common offsets; ordinary unnormalized weights can overflow mathematically.

## Untranslated computational functions

- `loo_moment_match`: full model-callback moment-matching or subsampling workflow remains unimplemented.
- `loo_moment_match.default`: full model-callback moment-matching or subsampling workflow remains unimplemented.
- `loo_subsample`: full model-callback moment-matching or subsampling workflow remains unimplemented.
- `loo_subsample.function`: full model-callback moment-matching or subsampling workflow remains unimplemented.

## Excluded explicit exports

- `.ndraws`: Draw-count metadata getter.
- `example_loglik_array`: Bundled example dataset accessor.
- `example_loglik_matrix`: Bundled example dataset accessor.
- `extract_log_lik`: Stan object extraction/integration without a numerical kernel.
- `find_model_names`: R model-name metadata extraction.
- `gpdfit`: Pass-through to posterior::gpdfit; reused from sibling posterior, not counted as a loo computation.
- `is.kfold`: R class predicate only.
- `is.loo`: R class predicate only.
- `is.psis`: R class predicate only.
- `is.psis_loo`: R class predicate only.
- `is.sis`: R class predicate only.
- `is.tis`: R class predicate only.
- `is.waic`: R class predicate only.
- `kfold`: Bare external-method dispatch generic: loo supplies no fitting/computational kfold method.
- `nlist`: R named-list construction.
- `obs_idx`: Subsample observation-index metadata getter.
- `pareto_k_influence_values`: Existing diagnostic vector getter.
- `pareto_k_values`: Existing diagnostic vector getter.
- `pointwise`: Existing pointwise matrix getter.
- `print_dims`: Presentation-only dimensions output.
- `psis_n_eff_values`: Existing diagnostic vector getter.
