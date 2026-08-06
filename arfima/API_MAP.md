# API map

| R `arfima` routine | Fortran routine or representation | Notes |
|---|---|---|
| `arfima` | `fit_arfima`, `fit_arfima_modes`, `arfima_estimate` | Typed model specification and result objects replace S3 lists. |
| `arfima0` | `arfima0_fit`, compatibility `arfima0` | Simplified exact MLE. |
| `arfima.sim` | `arfima_simulate`, compatibility `arfima_sim` | Optional supplied innovations and initial integrated values. |
| `sim_from_fitted` | `sim_from_fitted` | Static-regression fits supported; dynamic transfer simulation remains unsupported as upstream also documents. |
| `predict.arfima` | `arfima_forecast`, `predict_ARFIMA`, `predict_from_fitted` | Exact conditional covariance and integer-integration propagation. |
| `lARFIMA` | `dl_loglikelihood` plus `tacvf_arfima`; compatibility `lARFIMA` | Exact concentrated Gaussian likelihood. |
| `lARFIMAwTF` | compatibility `lARFIMAwTF` | Uses typed `transfer_spec`. |
| `tacvfARFIMA` | `tacvf_arfima`, compatibility `tacvfARFIMA` | Ordinary and seasonal long-memory mixtures. |
| `tacvfARMA` | `tacvf_arma`, compatibility `tacvfARMA` | Exact linear equations, arbitrary orders. |
| `tacvfFDWN` | `tacvf_fdwn`, compatibility `tacvfFDWN` | Gamma-function initial value and recurrence. |
| `tacvfFGN` | `tacvf_fgn`, compatibility `tacvfFGN` | Exact second-difference formula. |
| `tacvfHD` | `tacvf_pla`, compatibility `tacvfHD` | Borwein zeta approximation retained. |
| `tacvf` | `tacvf_from_fitted` | Returns a numeric autocovariance vector. |
| `ARToPacf` | `ar_to_pacf`, compatibility `ARToPacf` | Barndorff-Nielsen-Schou recursion. |
| `PacfToAR` | `pacf_to_ar`, compatibility `PacfToAR` | Stable inverse recursion. |
| `IdentInvertQ` | `identifiable_invertible`, compatibility `IdentInvertQ` | Stability/range checks plus numerical Fisher information. |
| `iARFIMA` | `arfima_information`, compatibility `iARFIMA` | Numerical spectral-information calculation. |
| `integ` | `integrate_series`, compatibility `integ` | Ordinary and seasonal inverse differencing. |
| `psiwts`, `wtsforexact` | `psi_weights`, `exact_integration_weights`; compatibility wrappers | Polynomial/filter weights. |
| transfer-function native `tfcalc` | `apply_transfer_function` | Multiple inputs and concatenated delta/omega vectors. |
| native `durlevsim` | `dl_simulate` | Arbitrary supplied standardized innovations. |
| `distance` | `mode_distance` | p-norm in natural parameter space. |
| `weed` | `weed_modes` | Removes nearby modes. |
| `bestModes` | `best_modes` | Sorts by likelihood and retains the requested count. |
| `removeMode` | `remove_mode` | Removes one mode. |
| `AIC`, `BIC`, `logLik`, `coef`, `vcov`, `residuals`, `fitted` | fields of `arfima_fit_result` | No R method dispatch required. |
| plotting and print methods | omitted | No computational kernel. |
| `arfimachanges` | omitted | Package-news text only. |
