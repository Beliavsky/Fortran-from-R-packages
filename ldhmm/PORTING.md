# Porting notes

## R-to-Fortran mapping

| R function | Fortran procedure |
|---|---|
| `ecld` | `ecld_create` |
| `ecld.pdf`, `ecld.cdf`, `ecld.ccdf` | `ecld_pdf`, `ecld_cdf`, `ecld_ccdf` |
| `ecld.mean`, `ecld.sd`, `ecld.var` | `ecld_mean`, `ecld_sd`, `ecld_variance` |
| `ecld.skewness`, `ecld.kurtosis` | `ecld_skewness`, `ecld_kurtosis` |
| `ldhmm` | `ldhmm_create` |
| `ldhmm.gamma_init` | `ldhmm_gamma_init` |
| `ldhmm.n2w`, `ldhmm.w2n` | `ldhmm_natural_to_working`, `ldhmm_working_to_natural` |
| `ldhmm.state_ld`, `ldhmm.state_pdf` | `ldhmm_state_distribution`, `ldhmm_state_pdf` |
| `ldhmm.mllk` | `ldhmm_mllk` |
| `ldhmm.log_forward`, `ldhmm.log_backward` | `ldhmm_log_forward`, `ldhmm_log_backward` |
| `ldhmm.decoding`, `ldhmm.viterbi` | `ldhmm_decode`, `ldhmm_viterbi` |
| `ldhmm.conditional_prob` | `ldhmm_conditional_prob` |
| Forecast functions | `ldhmm_forecast_*` |
| `ldhmm.pseudo_residuals` | `ldhmm_pseudo_residuals` |
| State-statistics functions | `ldhmm_ld_stats`, `ldhmm_calc_stats_from_obs`, `ldhmm_decode_stats_history` |
| Simulation functions | `ldhmm_simulate_state_transition`, `ldhmm_simulate_abs_acf` |
| `ldhmm.mle` | `ldhmm_fit` |
| `ldhmm.sma`, `ldhmm.ts_abs_acf` | `ldhmm_sma`, `ldhmm_abs_acf` |

## Numerical implementation changes

- Forward, backward, and Viterbi recursions are implemented directly in log
  space rather than by repeated rescaling. This preserves the same probability
  model while improving underflow resistance.
- Generalized-normal random values are generated using a gamma variate and a
  random sign. The R code used inverse-CDF root finding on a fixed interval.
- `ldhmm_create` calculates stationary `delta` when it is omitted. The R
  constructor can retain `NULL` until a later working-to-natural transform.
- The one-state likelihood uses the mathematically correct negative sum of log
  densities. Upstream 0.6.1 appears to use negative sum of densities in that
  special branch.
- Annualized allocation is recalculated for every state. Upstream 0.6.1 appears
  to update only the final row in that branch.
- Pseudo-residual integration uses a range expanded by 25 percent of the data
  range on each side, avoiding sign-dependent grid truncation.
- NaN observations are treated as missing observations in all likelihood
  recursions, with emission factor one.

## Optimization

The R package delegates to `nlm` or `optimx`. This translation supplies
self-contained finite-difference BFGS and Nelder-Mead implementations. The
optimizer name `nlm` maps to Nelder-Mead; `optimx` maps to BFGS.

## Excluded R infrastructure

Plotting, graphics callbacks, `xts`/`zoo` conversion and indexing, FRED
network access, YAML data configuration, RData loading, and package-specific
archive readers are excluded. Generic price-to-log-return conversion remains
available as `ldhmm_prices_to_log_returns`.
