# Computational API map

This map groups original R functions by numerical behavior. Plain Fortran procedures replace aliases, S3 dispatch, formulas, and report formatting.

## Direct or definition-preserving translations

| Original area | Fortran procedures |
|---|---|
| Return calculation and conversion | `calculate_returns`, `cumulative_return`, `annualized_return`, `annualized_excess_return`, `excess_returns`, `relative_returns`, `centered_returns`, `level_from_returns` |
| Portfolio accounting | `portfolio_returns`, `portfolio_result`, `wealth_index` |
| Descriptive and partial moments | procedures in `statistics_mod` |
| Drawdown calculations | procedures in `drawdown_mod` |
| Historical/Gaussian/modified VaR and ES | procedures in `risk_mod` |
| Lognormal VaR and ES | `lognormal_var`, `lognormal_es` |
| GPD tail fitting and risk | `gpd_fit`, `gpd_var_value`, `gpd_es_value` |
| Monte Carlo risk | `monte_carlo_asset_risk`, `monte_carlo_portfolio_risk` |
| Kernel portfolio risk | `kernel_portfolio_risk` |
| Moving-block risk errors | `bootstrap_risk_standard_errors` |
| CAPM/SFM and timing regressions | procedures in `capm_mod` |
| Rolling/expanding CAPM | `rolling_sfm`, `expanding_sfm` |
| Conditional/dynamic CAPM numerical core | `conditional_capm_fit` |
| Co-moment matrices and contractions | procedures in `comoments_mod` |
| M2/M3/M4 structured estimators | `structured_covariance`, `structured_coskewness`, `structured_cokurtosis`; independent M3 targets accept `unbiased_marg=.true.` |
| Exact M2/M3/M4 shrinkage | `exact_m2_shrinkage`, `exact_m3_shrinkage`, `exact_m4_shrinkage`, `exact_shrinkage_result` |
| Exact shrinkage loss kernels | `exact_vm2_terms`, `exact_vm3_terms`, `exact_vm3_kstat_terms`, `exact_vm4_terms`, `solve_shrinkage_qp` |
| Lightweight plug-in shrinkage | `shrink_covariance`, `shrink_coskewness`, `shrink_cokurtosis`, `multi_target_shrink_*` |
| EWMA co-moments | `ewma_covariance`, `ewma_coskewness`, `ewma_cokurtosis` |
| M3/M4 Moment Component Analysis | `m3_mca`, `m4_mca` |
| Nearest Comoment Estimation | `nearest_comoment_estimator`, `nce_result` |
| Portfolio Gaussian/modified moment risk | procedures in `portfolio_risk_mod` |
| Rolling and expanding statistics | procedures in `rolling_mod` |
| Capture and outperformance table values | `capture_ratios`, `outperformance_probabilities` |
| Common performance table values | `compute_performance_summary`, `performance_summary` |

## Numerical analogues and differences

| Original area | Fortran analogue |
|---|---|
| `Return.clean(method="boudt")` | Median/MAD winsorization in `clean_boudt` |
| `M2.shrink`, `M3.shrink`, `M4.shrink` | Exact analytical A/b construction and target families in `exact_m*_shrinkage`; projected-gradient simplex QP replaces `quadprog` |
| Constant-correlation M3/M4 targets | Standardized symmetry-class averaging |
| `MM.NCE` | Factor-model moment construction, PCA initialization, feasibility projection, and deterministic coordinate optimization |
| `NCEconstructW` | Identity, diagonal, ridge-diagonal, and ridge-identity order weighting rather than the complete full element weight matrix |
| `VaR.gpd`, `ES.gpd` intervals | Numerical Hessian and delta-method intervals rather than exact profile intervals |
| Monte Carlo and bootstrap | Reproducible xorshift/Box-Muller Fortran stream |
| Calendar rolling/rebalancing | Integer observation windows |

## Excluded infrastructure

- Plotting and interactive functions
- Formulas, S3/S4 methods, model frames, attributes, and `xts`/`zoo` metadata
- Formatted report/table rendering; numerical table values are exposed separately
- Exact optional external-package behavior
- Packaged datasets and R data import helpers

## Finite-sample correction notes

- `exact_vm3_kstat_terms` ports the original sixth-order k-statistic coefficients for the unbiased coskewness MSE path.
- The exact routines return the original diagnostic surface: target vectors, A, b, shrinkage intensities, sample estimate, and corrected estimate.
- Multiple observed factors add one one-factor target per factor, matching the R workflow.
- The all-target paths are covered by deterministic b-vector references in `test_finite_sample_corrections`.
