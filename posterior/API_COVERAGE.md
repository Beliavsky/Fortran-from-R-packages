# API coverage and export inventory

Status: **substantial**. **91 of 113 (80.53%)** exported computational R
functions have meaningful Fortran mappings. This fraction counts mapped
functions, including substantial and partial specializations; it is not a
measure of complete R compatibility. Every mapped function has one manifest
entry, with exact case and punctuation and public implementing procedures.

## Basis

The pinned NAMESPACE has 140 distinct explicit `export(...)` names.
We retain 113 computational names and exclude only the 27 interfaces
listed below. Numeric constructors/converters, draw reordering/subsetting,
array shape operations, higher-order transforms and random-variable arithmetic
remain in the denominator even when unimplemented. Registered S3 methods are
not separately exported by name in this NAMESPACE and are not double-counted.
`exportMethods(diag)` and `exportMethods(drop)` repeat already-counted exported
generics. Diagnostic aliases such as `summarize_draws` remain distinct exports.

## Mapped functions

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| %**% | posterior_draws | rvar_matmul | substantial |
| %in% | posterior_transforms | rvar_in | substantial |
| E | posterior_draws | rvar_expectation | substantial |
| Pr | posterior_draws | rvar_probability | substantial |
| as_draws_array | posterior_draws | matrix_to_draws | substantial |
| as_draws_matrix | posterior_draws | draws_to_matrix | substantial |
| as_rvar | posterior_draws | rvar_constant | partial |
| as_rvar_integer | posterior_transforms | rvar_as_integer | substantial |
| as_rvar_logical | posterior_transforms | rvar_as_logical | substantial |
| autocorrelation | posterior_diagnostics | autocorrelation | substantial |
| autocovariance | posterior_diagnostics | autocovariance | substantial |
| bind_draws | posterior_draws | bind_draws_arrays | substantial |
| cdf | posterior_draws | rvar_cdf | substantial |
| diag | posterior_transforms | rvar_diagonal, rvar_diagonal_matrix | substantial |
| dissent | posterior_predictive | dissent | substantial |
| draws_array | posterior_draws | matrix_to_draws | substantial |
| entropy | posterior_predictive | entropy | substantial |
| ess_basic | posterior_diagnostics | ess_basic | substantial |
| ess_bulk | posterior_diagnostics | ess_bulk | substantial |
| ess_mean | posterior_diagnostics | ess_mean | substantial |
| ess_median | posterior_diagnostics | ess_median | substantial |
| ess_quantile | posterior_diagnostics | ess_quantile | substantial |
| ess_sd | posterior_diagnostics | ess_sd | substantial |
| ess_tail | posterior_diagnostics | ess_tail | substantial |
| extract_variable | posterior_draws | extract_draws_variable | substantial |
| extract_variable_matrix | posterior_draws | extract_draws_variable_matrix | substantial |
| for_each_draw | posterior_transforms | for_each_draw | partial |
| gpdfit | posterior_pareto | gpdfit | substantial |
| is_constant | posterior_diagnostics | is_constant | partial |
| mad | posterior_draws | rvar_over_draws_summary | partial |
| match | posterior_transforms | rvar_match | substantial |
| mcse_mean | posterior_diagnostics | mcse_mean | substantial |
| mcse_median | posterior_diagnostics | mcse_median | substantial |
| mcse_quantile | posterior_diagnostics | mcse_quantile | substantial |
| mcse_sd | posterior_diagnostics | mcse_sd | substantial |
| merge_chains | posterior_draws | merge_draws_chains | substantial |
| modal_category | posterior_predictive | modal_category | substantial |
| mutate_variables | posterior_transforms | mutate_variables | partial |
| order_draws | posterior_draws | order_draws_array | substantial |
| pareto_convergence_rate | posterior_pareto | pareto_convergence_rate | substantial |
| pareto_diags | posterior_pareto | pareto_diags | substantial |
| pareto_khat | posterior_pareto | pareto_khat | substantial |
| pareto_khat_threshold | posterior_pareto | pareto_khat_threshold | substantial |
| pareto_min_ss | posterior_pareto | pareto_min_ss | substantial |
| pareto_pit | posterior_predictive | pareto_pit | substantial |
| pareto_smooth | posterior_pareto | pareto_smooth | substantial |
| pgeneralized_pareto | posterior_pareto | pgeneralized_pareto | substantial |
| pit | posterior_predictive | pit | substantial |
| ps_convergence_rate | posterior_pareto | ps_convergence_rate | substantial |
| ps_khat_threshold | posterior_pareto | ps_khat_threshold | substantial |
| ps_min_ss | posterior_pareto | ps_min_ss | substantial |
| ps_tail | posterior_pareto | ps_tail | substantial |
| ps_tail_length | posterior_pareto | ps_tail_length | substantial |
| qgeneralized_pareto | posterior_pareto | qgeneralized_pareto | substantial |
| quantile2 | posterior_diagnostics | quantile2 | substantial |
| r_scale | posterior_diagnostics | r_scale | substantial |
| resample_draws | posterior_draws | resample_draws_array, subset_draws_array | partial |
| rhat | posterior_diagnostics | rhat | substantial |
| rhat_basic | posterior_diagnostics | rhat_basic | substantial |
| rhat_nested | posterior_diagnostics | rhat_nested | substantial |
| rvar_all | posterior_draws | rvar_summary | substantial |
| rvar_any | posterior_draws | rvar_summary | substantial |
| rvar_apply | posterior_transforms | rvar_apply_margin | partial |
| rvar_ifelse | posterior_draws | rvar_ifelse | substantial |
| rvar_is_finite | posterior_draws | rvar_special_masks | substantial |
| rvar_is_infinite | posterior_draws | rvar_special_masks | substantial |
| rvar_is_na | posterior_draws | rvar_special_masks | partial |
| rvar_is_nan | posterior_draws | rvar_special_masks | partial |
| rvar_mad | posterior_draws | rvar_summary | substantial |
| rvar_max | posterior_draws | rvar_summary | substantial |
| rvar_mean | posterior_draws | rvar_summary | substantial |
| rvar_median | posterior_draws | rvar_summary | substantial |
| rvar_min | posterior_draws | rvar_summary | substantial |
| rvar_prod | posterior_draws | rvar_summary | substantial |
| rvar_quantile | posterior_draws | rvar_quantiles | substantial |
| rvar_range | posterior_draws | rvar_ranges | substantial |
| rvar_sd | posterior_draws | rvar_summary | substantial |
| rvar_sum | posterior_draws | rvar_summary | substantial |
| rvar_var | posterior_draws | rvar_summary | substantial |
| sd | posterior_draws | rvar_over_draws_summary | partial |
| split_chains | posterior_draws | split_draws_chains | substantial |
| subset_draws | posterior_draws | subset_draws_array | substantial |
| summarise_draws | posterior_diagnostics | summarise_draws | partial |
| summarize_draws | posterior_diagnostics | summarise_draws | partial |
| thin_draws | posterior_draws | thin_draws_array | partial |
| u_scale | posterior_diagnostics | u_scale | substantial |
| uniformity_test | posterior_predictive | uniformity_test | substantial |
| var | posterior_draws | rvar_over_draws_summary | partial |
| variance | posterior_draws | rvar_over_draws_summary | partial |
| weight_draws | posterior_draws | weight_draws_array, extract_draws_weights | partial |
| z_scale | posterior_diagnostics | z_scale | substantial |

The manifest notes define per-mapping scope and omitted behavior.

## Untranslated computational exports (22)

- `as_draws`
- `as_draws_df`
- `as_draws_list`
- `as_draws_rvars`
- `as_rvar_factor`
- `as_rvar_numeric`
- `as_rvar_ordered`
- `draws_df`
- `draws_list`
- `draws_matrix`
- `draws_rvars`
- `drop`
- `extract_list_of_variable_arrays`
- `extract_variable_array`
- `rdo`
- `repair_draws`
- `rfun`
- `rstar`
- `rvar`
- `rvar_factor`
- `rvar_ordered`
- `rvar_rng`

These include richer R-list/data-frame/rvar conversion and construction,
factor/ordered-level handling, general variable-array reconstruction, automatic
format dispatch, dimension dropping/repair, R closure/expression wrapping,
random-number-generator wrapping, and the classifier-based Rstar diagnostic.
They are not silently classified as presentation-only interfaces.

## Excluded interfaces, individually justified

| R export | Exclusion reason |
|---|---|
| chain_ids | Retrieves chain metadata IDs. |
| default_convergence_measures | Returns diagnostic function-name strings. |
| default_mcse_measures | Returns MCSE function-name strings. |
| default_summary_measures | Returns default summary function-name strings. |
| draw_ids | Retrieves draw metadata IDs. |
| draws_of | Retrieves the R rvar backing-draw slot. |
| draws_of<- | Assigns the R rvar backing-draw slot and metadata. |
| example_draws | Example dataset retrieval. |
| is_draws | R class/type predicate, with no numerical draw operation. |
| is_draws_array | R class/type predicate, with no numerical draw operation. |
| is_draws_df | R class/type predicate, with no numerical draw operation. |
| is_draws_list | R class/type predicate, with no numerical draw operation. |
| is_draws_matrix | R class/type predicate, with no numerical draw operation. |
| is_draws_rvars | R class/type predicate, with no numerical draw operation. |
| is_rvar | R class/type predicate, with no numerical draw operation. |
| is_rvar_factor | R class/type predicate, with no numerical draw operation. |
| is_rvar_ordered | R class/type predicate, with no numerical draw operation. |
| iteration_ids | Retrieves iteration metadata IDs. |
| nchains | R representation dimension/metadata accessor. |
| ndraws | R representation dimension/metadata accessor. |
| niterations | R representation dimension/metadata accessor. |
| nvariables | R representation dimension/metadata accessor. |
| rename_variables | Variable-name metadata operation. |
| reserved_variables | Returns reserved metadata-variable names. |
| set_variables | Variable-name metadata operation. |
| variables | Retrieves variable names. |
| variables<- | R variable-name assignment only. |
