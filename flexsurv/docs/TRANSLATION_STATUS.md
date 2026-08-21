# Translation status

This document separates computational coverage from R-language infrastructure.

## Complete or substantially complete numerical translations

| Upstream area | Fortran implementation | Status |
|---|---|---|
| Main parametric distributions | `flexsurv_distributions` | Density, CDF/survival, hazard/cumhaz, quantile, RNG, means/RMST |
| Generalized gamma/F | `flexsurv_distributions` | Current and original parameterizations |
| `flexsurvreg` likelihood | `flexsurv_fit` | Exact/right/left/interval censoring, delayed entry, right truncation, weights, background hazards |
| Ancillary-parameter regression | `flexsurv_fit` | Covariates on every built-in distribution parameter |
| Source-style starts/optimizers | `flexsurv_fit` + supplied `survival_aft` | `survreg`-style AFT starts, flexsurv moment starts, BFGS, Nelder-Mead, fallback |
| Custom distributions | `flexsurv_custom` | Callback likelihood under censoring/truncation |
| Royston-Parmar distributions/fits | spline modules | Hazard/odds/normal scales, log/identity time, QP initialization |
| `spline="splines2ns"` | `flexsurv_splines2ns` | Natural cubic B-spline basis, first derivative, linear boundary extrapolation, fitting integration |
| Spline ancillary interactions | `flexsurv_spline_interactions` | Independent regressions on gamma0, gamma1, ... with either spline basis |
| `flexsurvmix` | `flexsurv_mixture_full` | Known/unknown/partial event types, component covariates, multinomial-logit membership, EM/direct MLE |
| Louis information | `flexsurv_mixture_full` | Analytic multinomial membership derivatives plus component-specific numerical derivatives |
| `fmixmsm` | `flexsurv_fmixmsm` | Path/final probabilities, means/quantiles, parameter-bootstrap CIs |
| Markov multi-state prediction | `flexsurv_multistate` | Kolmogorov forward ODE via supplied deSolve translation |
| Semi-Markov simulation | `flexsurv_multistate_uncertainty` | Clock-reset competing transitions, predictable time-dependent covariates, `pmatrix.simfs`, LOS, CIs |
| Shared multi-state regression covariance | `flexsurv_shared_multistate` | One shared coefficient/covariance vector, cross-transition hazard covariance, joint bootstrap |
| Final-state summaries | `flexsurv_final_states` | `pfinal_fmsm`, `simfinal_fmsm`, conditional times and bootstrap CIs |
| `msfit`/AJ comparison core | `flexsurv_ajfit` plus multi-state modules | Parametric versus Nelson-Aalen/Aalen-Johansen helpers |
| Advanced `standsurv` | `flexsurv_standardize_advanced` | Rate tables, all-cause/attributable survival/hazard/RMST/quantile, transforms, contrasts, delta/bootstrap uncertainty |
| `survrtrunc` | `survrtrunc_fit` | Time-reversed estimator |
| `flexsurvrtrunc` | `fit_flexsurvrtrunc` | Joint/final likelihoods, fixed parameters, covariance/bootstrap |
| Fractional polynomials | `flexsurv_fracpoly` | Basis and derivative basis |
| Diagnostics | `flexsurv_diagnostics` | Cox-Snell residuals, hazard ratios, normal parameter draws |

## Residual numerical differences

The remaining numerical differences are now narrow rather than missing model classes:

1. The Fortran `tcovs` interface operates on explicit model-matrix column rates
   (and direct spline-gamma rates) instead of rebuilding an R formula/model matrix
   from raw covariates at every transition.  Linear predictable covariates and
   their interaction columns are representable exactly when their design-column
   rates are supplied.
2. Source-style `survreg` starts are used for right-censored/event-only AFT cases.
   R's `survreg` handles a wider set of response/formula edge cases than the
   supplied survival AFT translation, so unusual interval/stratified initialization
   can still start from a different point.  The fitted flexsurv likelihood itself
   is unchanged.
3. Component survival derivatives in Louis information are obtained through the
   supplied `numDeriv` implementation.  Membership derivatives/Hessians are
   analytic.  This mirrors the upstream split but floating-point step choices can
   still produce tiny covariance differences.
4. The native `splines2ns` implementation reproduces the numerical natural-spline
   construction used by flexsurv, but R attributes, column names and sparse-matrix
   representations are intentionally absent.

## Intentionally omitted non-computational infrastructure

- R formulas, model frames, environments and NSE
- S3 print/plot/lines/tidy/glance/augment dispatch
- ggplot2/tibble/dplyr/tidyr/purrr presentation plumbing
- R factors, dates, names/dimnames, recycling and warning semantics
- plotting and documentation/vignette build machinery
