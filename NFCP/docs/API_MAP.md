# R-to-Fortran API map

| R export | Fortran API | Notes |
|---|---|---|
| `NFCP_parameters` | `initialize_model`, `nfcp_parameterization_t%names` | Typed model plus R-order parameter-name generation |
| `NFCP_domains` | `default_parameter_bounds` | Same default parameter classes and practical bounds |
| `A_T` | `nfcp_a_t` | Elemental scalar function; array syntax works naturally |
| `cov_func` | `nfcp_covariance` | Dense N by N state covariance |
| `stitch_contracts` | `stitch_contract_numbers`, `stitch_by_maturity` | Numeric outputs and selected source-column indices |
| `NFCP_Kalman_filter` | `nfcp_kalman_filter` | Missing data use IEEE NaN |
| `NFCP_MLE` | `nfcp_fit_mle` | Differential evolution plus local L-BFGS-B, not `rgenoud` |
| `spot_price_forecast` | `spot_price_forecast` | Median plus optional percentile columns |
| `futures_price_forecast` | `futures_price_forecast` | Median plus optional percentile columns |
| `spot_price_simulate` | `spot_price_simulate` | Optional antithetic paths and deterministic seed |
| `futures_price_simulate` | `futures_price_simulate` | State, futures, and spot outputs |
| `TSfit_volatility` | `tsfit_volatility` | Theoretical and empirical annualized volatility |
| `European_option_value` | `european_option_value` | Value, volatility, and principal Greeks |
| `American_option_value` | `american_option_value` | Native Longstaff-Schwartz implementation |

## Representation changes

R named vectors are represented by `nfcp_model_t`. Dynamic lists are replaced
by typed result structures. Contract data and time-to-maturity data are always
passed as matrices with shape `(observation, contract)`. Missing observations
are IEEE NaNs.

The package data object `SS_oil` remains in the retained original R archive;
it is not converted into a Fortran object because R's `.rda` serialization is
not a portable numerical interchange format.
