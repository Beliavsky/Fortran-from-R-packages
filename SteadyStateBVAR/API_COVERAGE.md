# API coverage

Coverage is based on distinct exported computational R functions. Plotting, printing-only, formatting, dataset, download, interactive, and presentation-only exports are excluded. A mapped function may still have material R-interface or sampler differences.

| R function | Fortran procedure | Status | Main differences |
|---|---|---|---|
| `bvar` | `bvar_create` | substantial | Numeric container only; no S3/ts/name semantics. |
| `setup` | `setup_model` | substantial | Numerical setup and OLS core translated; R object/list semantics omitted. |
| `priors` | `set_priors` | substantial | Homoscedastic and complete RW/AR1 SV prior parameter sets translated; R named-list handling and warning text omitted. |
| `restrict_beta` | `restrict_beta` | substantial | Same `1e-7` near-zero restriction variance; no R warnings. |
| `fit` | `fit_model` | substantial | Homoscedastic/RW/AR1 model targets and predictive recursions translated with native Gibbs/Metropolis-within-Gibbs MCMC instead of Stan/NUTS; no Stan diagnostics, adaptive NUTS, chain object, or exact RNG parity. |
| `forecast` | `forecast_model` | substantial | Predictive mean/median, type-7 intervals, and annualized-growth transformation translated; no plotting/ts metadata/steady-state display overlay. |
| `conditional_forecast` | `conditional_forecast_model` | substantial | Homoscedastic equality-conditioned Gaussian forecasts translated with an equivalent projection formulation; upstream also rejects SV conditional forecasts. |
| `IRF` | `irf_model` | substantial | OIRF/GIRF, type-7 intervals, growth transformation, and time-specific SV covariance selection translated; response/impulse plot selection and plotting omitted. |
| `ppi` | `ppi` | complete | Numeric operation translated; R list return replaced by output arguments. |
| `summary.bvar` | `bvar_summary` | substantial | Mean/median beta/Psi/Sigma plus `A`, RW `phi`, and AR1 `gamma_0`/`gamma_1`/`Phi` translated, including SV `t` selection; R printing/filtering/rounding/labels omitted. |

Package translation status: **substantial**. Mapped computational functions: **10 of 10 (100%)**. That percentage measures function mapping coverage, not full R/Stan compatibility.

## Remaining compatibility differences

The principal remaining differences are execution-layer rather than missing model equations:

- the upstream fitter uses Stan/NUTS, while the Fortran fitter uses native fixed-tuning Gibbs/Metropolis-within-Gibbs updates;
- no Stan warmup/adaptation diagnostics, multiple-chain object, divergent-transition reporting, effective sample size, or R-hat metadata are produced;
- intrinsic Fortran RNG streams do not match Stan/R RNG streams;
- R S3 objects, names, `ts` attributes, printing, plotting, and data-frame/list wrappers are intentionally outside the computational translation;
- proposal-scale tuning for SV latent states and AR1 gamma parameters is caller-controlled rather than automatically adapted.
