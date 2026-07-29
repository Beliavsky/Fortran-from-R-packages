# API mapping

## Model fitting and simulation

| Original R entry point | Fortran API | Notes |
|---|---|---|
| `fit_ou_nonlinear_tmg` | `fit_ou_nonlinear_tmg` | Native single-level fit and diagnostics. |
| `fit_ou_nested` | `fit_ou_nested` | Levels 1, 2, and 3; level switches and value anchor. |
| `fit_ou_nested_mi` | `fit_ou_nested_mi` | Three-dimensional imputation array and Rubin pooling. |
| Stan generative equations | `simulate_ou_nested` | Single-, two-, and three-level simulation. |
| Stan pointwise likelihood | `ou_log_likelihood` | Total and optional pointwise log likelihood. |
| Stan mean equation | `ou_mean_increment` | One-step conditional mean increment. |

## Data preparation and configuration

| Original | Fortran |
|---|---|
| `zscore_train` | `zscore_train` |
| `compute_common_factor` | `compute_common_factor` |
| TMG orthogonalization | `orthogonalize_series` |
| weighted COM transformed data | `weighted_com_statistics` |
| `ou_level_spec` | `ou_level_spec` |
| column-name matching | `align_columns_indices` |

## Diagnostics and model comparison

| Original | Fortran |
|---|---|
| `evaluate_oos` | `evaluate_oos` |
| `evaluate_oos_nested` | `evaluate_oos_nested` |
| `compare_models_loo` | `compare_models_loo` |
| `count_divergences` | `count_divergences` |
| `validate_ou_fit` | `validate_ou_fit` |
| `extract_convergence_evidence` | `extract_convergence_evidence` |
| `kappa_stability_evidence` | `kappa_stability_evidence` |
| internal LOO computation | `psis_loo` |
| Rubin pooling | `rubin_combine` |

## Extractors and numerical plot data

| Original | Fortran |
|---|---|
| `extract_posterior_summary` | `extract_posterior_summary` |
| `build_beta_tmg_table` | `build_beta_tmg_table` |
| `summarize_sv_sigmas` | `summarize_sv_sigmas` |
| `drift_decomposition_grid` | `drift_decomposition_grid` |
| `build_accounting_block` | `build_accounting_block` |
| `extract_mu_trajectory` | `extract_mu_trajectory` |
| `plot_beta_tmg` numerical data | `build_beta_tmg_table` |
| `plot_drift_curves` numerical data | `drift_decomposition_grid` |
| `plot_sv_evolution` numerical data | `fit%h_median` and `summarize_sv_sigmas` |

## Geometry engine

| Original | Fortran |
|---|---|
| `ou_geom_target` | `ou_geom_target` |
| `ou_geom_metric_euclidean` | `ou_geom_metric_euclidean` |
| `ou_geom_metric_riemannian` | `ou_geom_metric_riemannian` |
| `ou_geom_hmc` | `ou_geom_hmc` |
| E-BFMI helper | `ou_geom_ebfmi` |
| metric evaluation | `ou_geom_mass` |

The Fortran target stores procedure pointers for log probability, gradient, and
optional Hessian. The SoftAbs metric uses the supplied Hessian or a finite-
difference Hessian of the gradient.

## Typed results

Important public types include:

- `ou_input`
- `ou_options`
- `ou_priors`
- `ou_summary`
- `ou_draws`
- `ou_fit_result`
- `ou_diagnostics`
- `loo_result`
- `oos_metric`
- `stability_result`
- `ou_mi_result`
- `ou_geom_target_type`
- `ou_geom_metric_type`
- `ou_geom_hmc_result`

## Presentation/integration functions

The following original functions are not separate Fortran numerical routines:

- `plot_beta_tmg`, `plot_drift_curves`, and `plot_sv_evolution`: replaced by
  numeric outputs listed above.
- `export_model_comparison`: ordinary file presentation; model comparison is
  available through `compare_models_loo`.
- `ou_geom_bridge`: Stan/cmdstan object plumbing; direct Fortran procedure
  pointers remove the need for a bridge.
- `ou_nested_stan_code`: the original Stan file is retained at
  `original/inst/stan/ou_nested.stan`.
