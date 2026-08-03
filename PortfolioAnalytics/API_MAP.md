# Computational API map

## Direct or close numerical counterparts

| Upstream R function/family | Fortran counterpart |
|---|---|
| `portfolio.spec`, box/weight-sum constraints | `portfolio_constraints`, `initialize_constraints` |
| group/factor/turnover/diversification/position/return constraints | fields of `portfolio_constraints`, `check_constraints`, `constraint_violation` |
| `equal.weight` | `equal_weight_portfolio` |
| `inverse.volatility.weight` | `inverse_volatility_portfolio` |
| `turnover`, `var.portfolio`, `HHI`, `diversification` | same-named numerical routines |
| `random_portfolios`, `rp_sample`, `rp_grid` | `random_portfolios`, `random_grid_portfolios` |
| `optimize.portfolio` | `optimize_portfolio` |
| mean-variance / maximum-return / quadratic-utility paths | objective constants plus projected-gradient solver |
| random / DE-style optimization | random-search and differential-evolution solvers |
| ES, STARR, risk-budget, semideviation, drawdown objectives | objective constants and `evaluate_portfolio_objective` |
| `extract_risk` ES/CSM/EQS calculations | `historical_es`, `conditional_second_moment`, `expected_quadratic_shortfall` |
| efficient-frontier functions | `create_efficient_frontier` |
| `optimize.portfolio.rebalancing` | `optimize_rebalancing` |
| `black.litterman` | `black_litterman` |
| `EntropyProg` | `entropy_pool` |
| `meucci.moments`, `meucci.ranking` | `meucci_moments`, `meucci_ranking` |
| `ac.ranking`, `centroid` | `ac_ranking`, `centroid` |
| `statistical.factor.model` | `fit_statistical_factor_model` |
| `covarianceSF/MF`, `coskewnessSF/MF`, `cokurtosisSF/MF` | corresponding Fortran routines |
| sample moment functions | `sample_moments`, `sample_coskewness`, `sample_cokurtosis` |
| robust moment helpers | `winsorize_returns`, `covariance_shrinkage`, `robust_covariance_huber` |

## Adapted interfaces

- The upstream package dispatches to ROI, CVXR, quadprog, GLPK, Symphony,
  DEoptim, GenSA, pso, mco, OSQP, and other solvers. The Fortran port instead
  supplies native projected-gradient, differential-evolution, random-search,
  and local-pattern algorithms.
- Position limits are handled with cardinality repair plus nonlinear search.
  This is useful for small and medium problems but is not an exact MILP solver.
- CSM and EQS are evaluated directly through scalar convex/tail calculations;
  no conic-programming backend is required.
- Rebalancing uses an explicit logical rebalance mask and integer rolling window
  rather than `xts` calendar endpoints.
- Objective and constraint specifications are typed fields rather than dynamic
  R lists and function names.

## Omitted R infrastructure

Plotting and chart functions, S3 print/summary/extract methods, `xts`/`zoo`
index handling, formula/data-frame interfaces, packaged datasets, foreach
parallelism, vignettes/reports, and dynamic R callbacks are not translated.
The exact external robust-covariance wrappers, regime object hierarchy,
multi-layer S3 portfolio composition, and solver-plugin objects are retained in
the upstream snapshot but are not represented one-for-one.
