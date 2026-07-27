# API Map

This file maps the computational surface of the original `MSGARCH` package to the modern Fortran project. R classes and display methods are intentionally omitted.

## Exported R operations

| Original operation | Fortran equivalent | Status |
|---|---|---|
| `CreateSpec` | `create_spec`, `regime_spec`, `msgarch_spec` | Implemented |
| `FitML` | `fit_ml` | Implemented numerical analogue |
| `FitMCMC` | `fit_mcmc` | Implemented numerical analogue |
| `DIC` | `dic_from_draws`, `mcmc_result%dic` | Implemented |
| `State` | `hamilton_filter` predicted/filtered/smoothed arrays and Viterbi path | Implemented |
| `TransMat` | `msgarch_spec%transition`, `transition_matrix_power`, `forecast_state_probabilities` | Implemented |
| `Volatility` | `conditional_volatility`, `posterior_volatility` | Implemented |
| `UncVol` | `unconditional_variances`, `unconditional_volatility`, `posterior_unconditional_volatility` | Implemented numerical analogue |
| `PredPdf` | `predictive_pdf`, `in_sample_pdf`, `predictive_distribution_forecast`, posterior predictive routines | Implemented |
| `PIT` | `pit_values`, `in_sample_cdf`, `predictive_cdf`, `predictive_distribution_forecast`, `posterior_pit` | Implemented |
| `Risk` | `risk_in_sample`, `risk_forecast`, `posterior_risk_forecast` | Implemented |
| `ExtractStateFit` | `extract_regime` plus ordinary fit/filter calls | Partial numerical equivalent |
| `simulate` | `simulate_msgarch`, `simulate_ahead` | Implemented |
| `predict` | `forecast_mean`, `forecast_volatility`, predictive distribution and state-probability routines | Implemented |
| `logLik` | `filter_result%loglik`, `fit_result%loglik` | Implemented |
| printing/summary/plotting | None | Excluded R infrastructure |

## Internal numerical components

| Original area | Fortran procedures | Status |
|---|---|---|
| Distribution densities | `innovation_pdf`, `innovation_logpdf` | Implemented for all six package distributions |
| Distribution CDF/quantile/RNG | `innovation_cdf`, `innovation_quantile`, `random_innovation` | Implemented |
| Distribution moments | `distribution_moments` | Implemented and used in stationarity checks |
| Model constraints | `regime_valid`, `spec_valid`, `parameter_bounds` | Implemented |
| Parameter transformations | `bounded_map`, `bounded_unmap`, `simplex_mapping`, `simplex_unmapping` | Implemented |
| Parameter fixing | `fit_ml(..., fixed_mask, fixed_values)` | Implemented |
| Regime-constant parameters | `fit_ml(..., tie_group)` | Implemented |
| Variance targeting | `variance_target_intercept` | Implemented using original formulas |
| Hamilton likelihood | `hamilton_filter` | Implemented |
| Smoothing | `hamilton_filter%smoothed` | Implemented |
| State decoding | `hamilton_filter%viterbi`, `viterbi_gaussian_hmm` | Implemented |
| Conditional volatility recursion | `initialize_state`, `update_state` | Implemented for all five models |
| ML optimization | `nelder_mead`, `fit_ml` | Numerical analogue |
| Inference | `numerical_hessian`, `invert_matrix` | Implemented numerical analogue |
| MCMC kernel | `log_posterior`, `fit_mcmc` | Numerical analogue |
| MCMC label sorting | `homogeneous_regimes`, `sort_parameters_by_variance` | Implemented |
| Posterior averaging | `msgarch_posterior` procedures | Implemented |
| Predictive KDE/empirical CDF | `predictive_distribution_forecast` | Implemented numerical analogue |
| Gaussian HMM initialization | `fit_gaussian_hmm` | Implemented numerical analogue |
| Gaussian mixture initialization | `fit_gaussian_mixture` | Implemented numerical analogue |

## Model names

The accepted model identifiers are case-insensitive after trimming:

- `sARCH`
- `sGARCH`
- `eGARCH`
- `gjrGARCH`
- `tGARCH`

## Distribution names

- `norm`
- `std`
- `ged`
- `snorm`
- `sstd`
- `sged`

## Not reproduced exactly

- The exact Vihola adaptive MCMC implementation
- Custom R sampler callbacks
- R's automatic starting-value search sequence and optimizer fallbacks
- R's class-preserving extraction methods
- R `density()` bandwidth selection and exact floating-point results
- R/Rcpp random-number streams
