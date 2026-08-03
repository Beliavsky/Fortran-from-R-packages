# PortfolioAnalytics-fortran

A self-contained modern Fortran/FPM translation of the reusable computational
core of the R package **PortfolioAnalytics 2.1.2**.

The API uses ordinary Fortran vectors, matrices, and derived types. Returns are
stored as `observations x assets`, matching the R package's conceptual layout.
Plotting, S3 objects, `xts`/`zoo` dates, data frames, report generation, and
parallel R orchestration are intentionally omitted.

## Main capabilities

- Portfolio constraints: box, weight sum, group, factor exposure, return target,
  leverage, turnover, diversification, transaction costs, and position counts.
- Objectives: minimum variance, maximum return, quadratic utility, maximum
  Sharpe, minimum historical ES, maximum STARR, risk parity, semideviation,
  maximum drawdown, HHI concentration, CSM, CSM ratio, and EQS.
- Solvers: projected gradient, differential evolution, random search, and local
  pattern refinement.
- Equal-weight, inverse-volatility, random, discretized-grid, and rolling
  rebalanced portfolios.
- Efficient-frontier construction and component volatility/ES calculations.
- Black-Litterman, entropy pooling, Meucci flexible-probability moments and
  ranking views, and Almgren-Chriss centroid ranking.
- PCA statistical factor models, covariance/coskewness/cokurtosis helpers,
  covariance shrinkage, winsorization, and Huber-style robust covariance.

## Build

```text
fpm build
fpm test
fpm run
```

On Windows with a recent GNU Fortran installation, the same commands work from
`cmd.exe` or PowerShell. The package has no external numerical-library or
runtime-backend requirement.

## Minimal example

```fortran
use portfolio_analytics

real(dp) :: returns(100,4)
type(portfolio_constraints) :: constraints
type(portfolio_options) :: options
type(portfolio_result) :: result

call initialize_constraints(constraints,4, &
  max_weight=[0.5_dp,0.5_dp,0.5_dp,0.5_dp])
options%objective = obj_min_variance
call optimize_portfolio(returns,constraints,options,result)
```

See `example/`, `app/demo_portfolio_analytics.f90`, `API_MAP.md`, and
`PORTING_NOTES.md` for complete details.
