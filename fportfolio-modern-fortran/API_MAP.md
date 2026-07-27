# Computational API map

## Direct or definition-preserving translations

| fPortfolio numerical surface | Modern Fortran API |
|---|---|
| `covEstimator` | `sample_estimator` |
| `lpmEstimator` | `lpm_estimator` |
| `slpmEstimator` | `slpm_estimator` |
| `kendallEstimator` | `kendall_estimator` |
| `spearmanEstimator` | `spearman_estimator` |
| covariance shrinkage | `shrinkage_estimator` |
| covariance portfolio risk | `covariance_risk` |
| historical VaR/CVaR | `historical_var`, `historical_es` |
| Normal VaR/ES | `normal_var`, `normal_es` |
| modified VaR/ES | `modified_var`, `modified_es` |
| covariance risk budgets | `covariance_risk_contributions` |
| Normal VaR/ES budgets | `normal_var_contributions`, `normal_es_contributions` |
| modified VaR/ES budgets | `modified_var_contributions`, `modified_es_contributions` |
| historical VaR/ES budgets | `historical_var_contributions`, `historical_es_contributions` |
| `pfolioReturn` | `portfolio_returns`, `geometric_portfolio_returns` |
| `pfolioMaxLoss` | `portfolio_max_loss` |
| drawdown risk | `drawdown_series`, `maximum_drawdown`, `drawdown_at_risk`, `conditional_drawdown_at_risk` |
| tail dependence matrix | `empirical_tail_dependence`, `normal_margin_tail_dependence` |
| `.cfgFit` | `cfg_tail_dependence` |
| feasible portfolio | `feasible_portfolio`, `project_feasible` |
| minimum variance | `minvariance_portfolio` |
| efficient portfolio | `efficient_portfolio` |
| maximum return | `maxreturn_portfolio` |
| tangency/max-ratio | `tangency_portfolio`, `maxratio_portfolio` |
| efficient frontier | `portfolio_frontier` |
| MAD portfolio | `minimum_mad_portfolio` |
| CVaR portfolio | `minimum_cvar_portfolio` |
| risk budgeting | `risk_parity_portfolio` |
| rolling portfolio workflows | `run_backtest` |
| weight smoothing | `exponential_weight_smoother` |
| net performance | `net_performance` |
| rolling risk statistics | `rolling_sigma`, `rolling_var`, `rolling_cvar`, `rolling_dar`, `rolling_cdar` |
| EMA/MACD/drawdown indicators | `ema_indicator`, `macd_indicator`, `drawdown_indicator` |
| rebalancing statistics | `rebalancing_statistics` |

## Tested numerical extensions

- Maximum-diversification optimization
- Exact subset cardinality and buy-in search for up to 22 assets
- Generic dense box/equality/inequality constrained projected QP and LP APIs
- Explicit fixed and proportional transaction costs
- Rolling stability distances
- CSV readers and command-line workflows

## Numerical analogues rather than backend reproductions

| Original path | Translation |
|---|---|
| `quadprog`, `ipop`, CLA, short-exact wrappers | projected dense convex QP with active constraint projections |
| GLPK MAD | smooth projected MAD optimization |
| GLPK CVaR | primal-dual finite-scenario CVaR optimization |
| `solnp`/`nlminb2` ratio optimization | projected analytical-gradient tangency optimization |
| external risk-budget optimizer | multiplicative projected risk-budget iterations |

## Explicit exclusions

- S4 classes, specifications, getters/setters, formulas, and `timeSeries` data
- Plotting, sliders, ternary maps, weight wheels, frontier plots, and reports
- AMPL, NEOS, GLPK, Symphony, Kestrel, `quadprog`, `ipop`, `Rsolnp`, and SOCP backend interfaces
- Exact generic mixed-integer LP or arbitrary nonlinear programming backend parity
- Imported MVE, MCD, OGK, NNVE, ARW, Bayesian-Stein, RMT, and external robust covariance implementations
- Imported generalized-spherical-Normal, GLD, and GH-tail dependency fits
- `bcp`, GARCH, wavelet, and multivariate-outlier monitoring wrappers
- Calendar/date class operations and exact R backtest object orchestration
- Exact R random streams, optimizer trajectories, or solver certificates
