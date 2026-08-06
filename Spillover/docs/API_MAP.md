# API map

| R function or dependency operation | Fortran procedure | Notes |
|---|---|---|
| `vars::VAR()` | `fit_var` | OLS VAR(p); supports none, constant, trend, and constant plus trend. Seasonal and exogenous regressors are omitted. |
| `vars::Phi()` | `ma_coefficients` | Returns `Phi(:,:,1)=I` through the requested horizon. |
| `g.fevd()` | `generalized_fevd`, `g_fevd` | Returns a `k x k x horizon` array rather than an R list of matrices. |
| `G.spillover()` | `generalized_spillover`, `g_spillover` | Returns `type(spillover_result)`. |
| `vars::fevd()` | `orthogonalized_fevd` | Cholesky FEVD, optionally for a supplied permutation. |
| `fastSOM::sot_avg_est()` | `orthogonal_average_sample` | Seeded Fisher-Yates permutation sampling. |
| `fastSOM::sot_avg_exact()` | `orthogonal_average_exact` | Direct lexicographic enumeration; guarded by `exact_limit`. |
| `O.spillover()` | `orthogonalized_spillover`, `o_spillover` | Supports `ortho_single`, `ortho_partial`, and `ortho_total`. |
| spillover data frame | `compatibility_table` | Produces the numeric `(k+2) x (k+1)` table layout. |
| `net()` | fields `to`, `from`, `net` in `spillover_result` | `net = to - from`. |
| `dynamic.spillover()` | `dynamic_spillover` | Numeric rolling output without dates or data frames. |
| `roll.spillover()` | `rolling_total_spillover`, `roll_spillover` | Generalized or orthogonalized rolling total index. |
| `total.dynamic.spillover()` | `total_dynamic_spillover` | Compatibility alias of rolling total connectedness. |
| `roll.net()` | `rolling_net_spillover`, `roll_net` | Rolling net connectedness matrix. |
| `plotdy()` | omitted | Plotting-only code. |
