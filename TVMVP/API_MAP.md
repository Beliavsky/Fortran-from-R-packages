# API map

| Upstream R routine | Fortran procedure/type | Status |
|---|---|---|
| `silverman` | `silverman` | Direct formula |
| `epanechnikov_kernel` | `epanechnikov_kernel` | Direct formula |
| `boundary_kernel` | `boundary_kernel` | Source-compatible and corrected modes |
| `two_fold_convolution_kernel` | `two_fold_convolution_kernel` | Exact Epanechnikov polynomial; numerical custom-kernel path |
| `local_pca` | `local_pca` | Algebraically equivalent smaller-Gram SVD |
| `localPCA` | `local_pca_all`, `localPCA` | Complete time-index result |
| `determine_factors` | `determine_factors` | Direct IC and reconstruction workflow |
| `residuals` | `factor_residuals` | Direct computation |
| `adaptive_poet_rho` | `adaptive_poet_rho` | Direct split-validation workflow |
| `estimate_residual_cov_poet_local` | same name | Direct POET/shrinkage workflow |
| `time_varying_cov` | `time_varying_cov` | Typed functional API |
| `sqrt_matrix` | `matrix_sqrt_abs` | Direct eigenvalue-absolute-value behavior |
| `compute_sigma_0` | `compute_sigma_0` | Direct tapered covariance |
| `compute_M_hat` | `compute_m_hat` | Direct statistic component |
| `compute_B_pT` | `compute_b_pt` | Direct statistic component |
| `compute_V_pT` | `compute_v_pt` | Direct statistic component |
| `hyptest` | `hyptest`, `hypothesis_result` | Bootstrap test |
| `comp_expected_returns` | `comp_expected_returns` | Same candidate orders; CSS rather than exact R ARIMA ML |
| `predict_portfolio` | `predict_portfolio`, `portfolio_prediction_result` | All three analytical portfolios |
| `expanding_tvmvp` | `expanding_tvmvp`, `expanding_window_result` | TV-MVP and equal-weight evaluation |
| `TVMVP`, `PortfolioPredictions`, `ExpandingWindow` | typed Fortran result structures | R6 state/printing replaced |
| plotting methods | omitted | Non-computational |
| `get_object_size` | omitted | R object-runtime utility |
