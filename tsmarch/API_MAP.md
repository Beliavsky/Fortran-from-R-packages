# API map

This document maps the main upstream computational concepts to the Fortran API.
R S3 accessors such as `coef`, `fitted`, `residuals`, `tscov`, and `tscor` map
to fields of the typed results rather than separate generic procedures.

## DCC

| Upstream concept | Fortran API |
|---|---|
| DCC model specification | `type(dcc_spec)` |
| DCC coefficients | `type(dcc_parameters)` |
| DCC filtering | `dcc_filter` |
| Constant-correlation filter | `dcc_constant_filter` |
| Marginal GARCH fits | `fit_marginal_garch` |
| DCC estimation | `estimate_dcc` |
| DCC simulation | `simulate_dcc_innovations` |
| DCC prediction | `forecast_dcc` |
| Stationarity measure | `dcc_stationarity` |
| Parameter vector conversion | `pack_dcc_parameters`, `unpack_dcc_parameters` |
| Conditional covariance/correlation | `dcc_fit%filtered%covariance`, `%correlation` |
| Likelihood/AIC/BIC | `dcc_fit%log_likelihood`, `%aic`, `%bic` |
| Numerical inference | `dcc_fit%hessian`, `%covariance`, `%standard_errors`, `%scores` |

## Copula-GARCH

| Upstream concept | Fortran API |
|---|---|
| Copula specification | `type(copula_spec)` |
| Marginal probability integral transform | `probability_transform` |
| Gaussian/Student copula transform | `copula_transform` |
| Copula filter | `copula_filter` |
| Copula estimation | `estimate_copula` |
| Copula prediction | `forecast_copula` |
| Uniform observations | `copula_fit%uniforms` |

## ICA and GO-GARCH

| Upstream concept | Fortran API |
|---|---|
| Data whitening | `whiten_data` |
| FastICA | `fastica` |
| RADICAL | `radical` |
| GO-GARCH specification | `type(gogarch_spec)` |
| GO-GARCH estimation | `estimate_gogarch` |
| GO-GARCH prediction | `forecast_gogarch` |
| Conditional covariance | `gogarch_covariance` or `gogarch_fit%covariance` |
| Conditional correlation | `gogarch_correlation` or `gogarch_fit%correlation` |
| Coskewness | `gogarch_coskewness` |
| Cokurtosis | `gogarch_cokurtosis` |
| Portfolio variance/skewness/kurtosis | `portfolio_variance`, `portfolio_skewness`, `portfolio_kurtosis` |

## Risk and distribution aggregation

| Upstream concept | Fortran API |
|---|---|
| Portfolio simulation paths | `portfolio_paths` |
| Simulation VaR/ES | `simulation_risk`, `value_at_risk`, `expected_shortfall` |
| Gaussian VaR/ES | `gaussian_risk` |
| Density convolution | `discrete_convolution`, `fft_convolution` |
| Density/CDF/quantile evaluation | `dfft`, `pfft`, `qfft` |

## Utilities and diagnostics

| Upstream concept | Fortran API |
|---|---|
| EWMA covariance | `ewma_covariance` |
| Ledoit-Wolf covariance | `lw_covariance` |
| Correlation-to-covariance | `cor2cov` |
| PSD repair | `make_psd` |
| Multivariate Gaussian/Student density | `multivariate_normal_density`, `multivariate_student_density` |
| Multivariate Gaussian/Student random draws | `rmvnorm`, `rmvt` |
| Fast combinations | `combn_fast` |
| Engle-Sheppard diagnostic | `escc_test` |
| Aggregate portfolio mean/volatility | `aggregate_mean`, `aggregate_sigma` |
| Lagged design matrices | `lag_matrix` |
| Triangular vector conversion | `lower_triangle`, `triangle_to_symmetric` |

## Status handling

Every major result has a `status` and `message` field.  Public status constants
are `tsm_success`, `tsm_invalid_argument`, `tsm_no_convergence`,
`tsm_numerical_failure`, and `tsm_singular`.
