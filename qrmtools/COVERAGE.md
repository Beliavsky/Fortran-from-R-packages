# Computational coverage

## Distribution and EVT routines

| Upstream area | Fortran implementation |
|---|---|
| `dGEV`, `pGEV`, `qGEV`, `rGEV` | `dgev`, `pgev`, `qgev`, `rgev` |
| `dGPD`, `pGPD`, `qGPD`, `rGPD` | `dgpd`, `pgpd`, `qgpd`, `rgpd` |
| `dPar`, `pPar`, `qPar`, `rPar` | `dpar`, `ppar`, `qpar`, `rpar` |
| GPD-tail functions | `dgpdtail`, `pgpdtail`, `qgpdtail`, `rgpdtail` |
| GEV quantile/PWM/MLE fitting | `fit_gev_quantile`, `fit_gev_pwm`, `fit_gev_mle` |
| GPD MOM/PWM/MLE fitting | `fit_gpd_mom`, `fit_gpd_pwm`, `fit_gpd_mle` |
| Hill estimator | `hill_estimator` |
| Mean-excess and tail functions | `mean_excess_np`, `mean_excess_gpd`, `tail_estimator_gpd` |
| Composite laws | `distribution_component` and `composite_*` |

Plotting wrappers around fitted shapes, tails, Hill estimates, and mean excess
were excluded; their numerical arrays and estimators are available.

## Risk measures and dependence bounds

| Upstream area | Fortran implementation |
|---|---|
| Nonparametric VaR/ES/RVaR | `var_np`, `es_np`, `rvar_np` |
| Student-t VaR/ES | `var_t`, `es_t`, `var_t01`, `es_t01` |
| GPD/Pareto/tail VaR and ES | `var_*`, `es_*` family |
| `gVaR`, `gEX` | `geometric_var`, `geometric_expectile` |
| `rearrange` | `rearrange_matrix` |
| `block_rearrange` | `block_rearrange_matrix` |
| `RA` | `ra_bounds` with a quantile callback |
| `ARA`, `ABRA` | `adaptive_ra_bounds`, with `use_blocks` selecting ABRA behavior |
| Crude homogeneous bounds | `crude_var_bounds_hom` |
| Wang/Pareto homogeneous bounds | `pareto_var_bounds_hom`, `wang_pareto_bounds` |
| Dual bound | `dual_bound_value`, `dual_worst_var` |

The Fortran adaptive interface unifies the upstream ARA and ABRA wrappers and
returns typed lower/upper rearrangement objects.

## Processes, pricing, diagnostics, and utilities

| Upstream area | Fortran implementation |
|---|---|
| `rBrownian`, `deBrowning` | `r_brownian`, `de_browning` |
| Black-Scholes and Greeks | `black_scholes`, `black_scholes_greeks` |
| `returns` | `compute_returns`, `invert_returns` |
| `hierarchical_matrix` | `hierarchical_matrix` and `hierarchy_node` |
| `alloc_ellip`, `conditioning`, `alloc_np` | corresponding allocation routines |
| `maha2_test`, `mardia_test` | corresponding diagnostic routines |
| `logLik_GARCH_11`, `fit_GARCH_11` | `loglik_garch_11`, `fit_garch_11` |
| `tail_index_GARCH_11` | `tail_index_garch_11` |

## Excluded infrastructure

- All plotting and graphical-device operations.
- `get_data` and network/data-provider adapters.
- R S3/list/data-frame formatting.
- `catch`, which is specific to R conditions.
- `fit_ARMA_GARCH`, which delegates to external `rugarch` and
  `HoltWinters` objects rather than implementing an independent algorithm.
- Compiled use of bundled R data objects.
