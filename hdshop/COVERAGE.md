# Computational coverage

## Exported HDShOP routine mapping

| HDShOP routine | Fortran implementation |
| --- | --- |
| `Sigma_sample_estimator` | `sigma_sample_estimator` |
| `CovShrinkBGP14` | `cov_shrink_bgp14`, `covshrinkbgp14` |
| `nonlin_shrinkLW` | `nonlin_shrink_lw`, `nonlin_shrinklw` |
| `InvCovShrinkBGP16` | `inv_cov_shrink_bgp16`, `invcovshrinkbgp16` |
| `mean_bs` | `mean_bs` |
| `mean_js` | `mean_js` |
| `mean_bop19` | `mean_bop19` |
| `MeanEstim` | `meanestim` |
| `CovarEstim` | `covarestim` |
| `new_MeanVar_portfolio` | `new_meanvar_portfolio` |
| `MeanVar_portfolio` | `meanvar_portfolio` |
| `validate_MeanVar_portfolio` | `validate_meanvar_portfolio` |
| `new_MV_portfolio_traditional` | `new_mv_portfolio_traditional` |
| `new_MV_portfolio_traditional_pgn` | `new_mv_portfolio_traditional_pgn` |
| `new_MV_portfolio_weights_BDOPS21` | `new_mv_portfolio_weights_bdops21` |
| `new_MV_portfolio_weights_BDOPS21_pgn` | `new_mv_portfolio_weights_bdops21_pgn` |
| `new_GMV_portfolio_weights_BDPS19` | `new_gmv_portfolio_weights_bdps19` |
| `new_GMV_portfolio_weights_BDPS19_pgn` | `new_gmv_portfolio_weights_bdps19_pgn` |
| `MVShrinkPortfolio` | `mvshrinkportfolio`, `mv_shrink_portfolio` |
| `test_MVSP` | `test_mvsp` |
| `plot_frontier` | `plot_frontier`, `bayesian_frontier` |
| `RandCovMtrx` | `randcovmtrx`, `random_covariance_matrix` |

The internal BOP19, BDOPS21, covariance-projection, alpha-variance, and Omega
formulas are also exposed in `hdshop_formulas` and `hdshop_inference`.

## Excluded presentation and data infrastructure

- S3 print, summary, and graphics methods.
- `ggplot2`, `lattice`, and bar-chart rendering.
- R data-frame/list dispatch and formula infrastructure.
- Direct loading of `SP_daily_asset_returns.RData`.

The numerical data for the efficient frontier and portfolio summaries are
returned in typed structures so callers can render them using their preferred
plotting library.
