# Upstream API map

This file maps the public R interface of MTS 1.2.1 to the Fortran API. A
"direct" entry means the numerical operation is represented directly, not that
R argument names, list layouts, printing, or interactive behavior are copied.

## Linear multivariate models

| R function | Fortran routine | Status |
|---|---|---|
| `VAR`, `VARs`, `refVAR` | `fit_var`, `fit_sparse_var`, `refine_var` | Direct matrix API |
| `VARorder`, `VARorderI` | `select_var_order`, `select_var_order_increasing` | Direct |
| `VARpred` | `predict_var` | Direct |
| `VARpsi` | `var_psi_weights` | Direct |
| `FEVdec` | `forecast_error_variance_decomposition` | Direct |
| `GrangerTest` | `granger_causality_test` | Direct |
| `VARMA`, `VMA`, `VMAe`, `refVARMA`, `refVMA`, `refVMAe` | `fit_varma`, `fit_vma` | Conditional-Gaussian matrix API; refinement is achieved through explicit re-fitting |
| `VARMACpp`, `VMACpp` | `varma_residuals`, `fit_varma`, `fit_vma` | Native Fortran replacement |
| `VARMAsim` | `simulate_varma` | Direct |
| `VARMApred` | `predict_varma` | Direct |
| `VARMAirf` | `varma_impulse_response` | Direct |
| `VARMAcov` | `varma_covariance` | Direct truncated MA representation |
| `PSIwgt` | `varma_psi_weights` | Direct |
| `PIwgt` | `pi_weight_matrices` | Direct polynomial recursion |
| `sVARMA`, `refsVARMA`, `sVARMApred` | `seasonal_lag_matrices` plus `fit_varma`/`predict_varma` | Seasonal polynomial composition is direct; specialized R sparse-parameter UI is not reproduced |
| `VARX`, `refVARX` | `fit_varx` | Direct matrix API |
| `VARXorder` | `select_varx_order` | Direct |
| `VARXpred` | `predict_varx` | Direct |
| `VARXirf` | `var_psi_weights` with fitted VARX autoregressive matrices | Computational equivalent |
| `Mlm` | `multivariate_linear_model` | Direct |
| `RLS` | `recursive_least_squares` | Direct |
| `REGts`, `refREGts`, `REGtspred` | `fit_regression_with_ar_errors` and explicit forecasts | Core estimation translated; R formula/list wrappers omitted |

## Diagnostics and transformations

| R function | Fortran routine | Status |
|---|---|---|
| `ccm` | `cross_correlation_matrices` | Direct |
| `mq`, `MTSdiag` | `multivariate_portmanteau`, `cross_correlation_matrices` | Core tests direct |
| `MarchTest`, `archTest` | `multivariate_arch_test`, `rank_arch_test` | Direct numerical tests |
| `MCHdiag` | `volatility_diagnostics` | Direct matrix API |
| `diffM` | `difference_matrix` | Direct |
| `msqrt` | `matrix_sqrt_symmetric` | Direct |
| `Vech`, `VechM` | `vech_matrix`, `unvech_matrix` | Direct |
| `Mtxprod`, `Mtxprod1` | `matrix_polynomial_product`, `matrix_polynomial_product_seasonal` | Direct generalized API |
| `Corner` | `corner_table` | Numerical Corner-style residual-scale table; no interactive display |
| `Eccm` | `extended_cross_correlation` | Numerical prefiltered ECCM table |
| `Kronid` | `approximate_kronecker_indices` | Approximate canonical-correlation identification |
| `Kronspec` | `kronecker_specification` | Direct 0/1/2 indicator construction |

## Volatility

| R function | Fortran routine | Status |
|---|---|---|
| `EWMAvol` | `ewma_covariance`, `fit_ewma_lambda` | Direct |
| `dccFit`, `dccPre` | `fit_dcc`, `dcc_correlations`, `dcc_log_likelihood` | Direct Engle/Tse-Tsui correlation layer; univariate margins are supplied by caller |
| `BEKK11` | `fit_bekk11`, `bekk11_filter`, `bekk11_log_likelihood` | Direct BEKK(1,1) core |
| `MCholV` | `modified_cholesky_volatility` | Direct computational counterpart |
| `comVol` | `common_volatility_components` | Direct eigen-analysis counterpart |
| `SCCor` | `constrained_group_correlation` | Direct grouped correlation constraints |
| `mtCopula` | multivariate Student-t density and RNG primitives | Full grouped t-copula optimizer is not implemented |

## Factors, Bayesian VAR, cointegration, and missing data

| R function | Fortran routine | Status |
|---|---|---|
| `hfactor` | `principal_components`, `constrained_factor_model` | Core factor estimators |
| `apca` | `asymptotic_pca` | Direct |
| `SWfore` | `stock_watson_forecast` | Direct matrix API |
| `BVAR` | `fit_bvar` | Conjugate matrix-form Bayesian VAR |
| `ECMvar`, `refECMvar` | `fit_vecm_known_beta` | Direct known-beta VECM core |
| `ECMvar1`, `refECMvar1` | `fit_vecm_johansen` | Johansen-style reduced-rank estimator |
| `Vmiss` | `estimate_missing_observation` | Direct numerical counterpart |
| `Vpmiss` | `estimate_partial_missing` | Direct numerical counterpart |

## Specialized upstream interfaces not reproduced one-for-one

The following routines encode large interactive or highly specialized model
specification systems. Their lower-level ingredients are present, but their R
list layouts and complete search/constraint logic are not direct ports:

- `Kronfit`, `refKronfit`, `Kronpred`
- `SCMid`, `SCMid2`, `SCMmod`, `SCMfit`, `refSCMfit`
- `tfm`, `tfm1`, `tfm2`, `Btfm2`
- scalar-ARIMA `backtest`
- grouped multivariate t-copula optimization in `mtCopula`

Plotting (`MTSplot` and plot branches), R data objects, console menus, formula
parsing, S3/list presentation methods, and Rcpp registration are outside the
Fortran scope.
