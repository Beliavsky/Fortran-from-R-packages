# R-to-Fortran API map

| R `wqc` routine | Modern Fortran routine/type | Status |
|---|---|---|
| `quantile_correlation_analysis(x, y, quantiles, wf, J, n_sim)` | `quantile_correlation_analysis(x, y, quantiles, wf, j_levels, n_sim, seed, stat, errmsg)` | Implemented |
| `apply_quantile_correlation(data, quantiles, wf, J, n_sim)` | `apply_quantile_correlation(data, quantiles, wf, j_levels, n_sim, seed, series_names, stat, errmsg)` | Implemented |
| Internal `QCSIS::qc` calls | `quantile_correlation`, delegating to `qcsis_mod::qc` | Implemented |
| `waveslim::mra(..., method="modwt", boundary="periodic")` | `waveslim::mra(..., method='modwt', boundary='periodic')` | Implemented through dependency |
| R `data.frame` pair result | `type(wqc_pair_result)` with `(level, quantile)` arrays | Implemented |
| Combined R `data.frame` | `type(wqc_multi_result)` containing one pair result per target | Implemented |
| `plot_quantile_heatmap` | None | Intentionally omitted |

## Result ordering

For a pair result, the first array dimension is wavelet level and the second is
quantile. Iterating levels outside and quantiles inside reproduces the row
ordering of the R data frame.
