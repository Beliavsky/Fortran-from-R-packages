# fPortfolio Modern Fortran

A modern Fortran translation of the computational portfolio-selection and
risk-analysis core of the R package `fPortfolio`.

The project uses plain arrays and derived result types. It intentionally does
not reproduce R S4 classes, formulas, plotting, `timeSeries` metadata, or
remote/external solver orchestration.

## License

The original package declares `GPL (>= 2)`. This translation uses
`SPDX-License-Identifier: GPL-2.0-or-later` in every Fortran source file.
`LICENSE` contains the complete GNU GPL version 2 text.

## Implemented numerical features

### Estimation

- Sample means and covariance
- Exponentially weighted means and covariance
- Lower-partial-moment and symmetric lower-partial-moment matrices
- Spearman and Kendall rank covariance estimators
- Diagonal-target covariance shrinkage
- Positive-definite matrix repair and standard dense linear algebra

### Portfolio risk

- Portfolio arithmetic and geometric return paths
- Covariance risk
- Historical VaR and Expected Shortfall
- Normal VaR and Expected Shortfall
- Cornish-Fisher modified VaR and modified Expected Shortfall
- Maximum loss
- Drawdown, maximum drawdown, DaR, and CDaR
- Covariance, Normal VaR, Normal ES, historical VaR/ES, and modified VaR/ES
  marginal/component risk and budgets
- Diversification ratio
- Empirical lower/upper tail dependence
- CFG tail-dependence coefficient after Normal marginal transformation

### Optimization

- Dense box- and linear-constrained quadratic programming
- Dense projected linear objectives
- Minimum-variance portfolios
- Target-return efficient portfolios
- Maximum-return portfolios
- Tangency/maximum-Sharpe portfolios
- Maximum-diversification portfolios
- Equal- or specified-budget risk-parity portfolios
- Mean absolute deviation portfolios
- Historical CVaR portfolios
- Efficient-frontier generation
- Exact small-universe cardinality/buy-in subset search
- Large-universe cardinality heuristic
- Budget, target-return, box, equality, and inequality constraints

### Backtesting and monitoring

- Rolling/fixed-window portfolio re-estimation
- Equal-weight, minimum-variance, tangency, risk-parity,
  maximum-diversification, and CVaR strategies
- Weight smoothing
- Proportional and fixed transaction costs
- Turnover, wealth, annualized return/volatility, Sharpe ratio, and drawdown
- Rolling sigma, VaR, CVaR, DaR, and CDaR
- EMA, MACD, drawdown indicators, turning points, rebalancing statistics,
  and rolling mean/covariance stability distances

## Build and test

Requirements:

- GNU Fortran or another Fortran 2018 compiler
- LAPACK
- BLAS

On Unix-like systems:

```sh
make debug
make release
```

The debug build enables bounds and runtime checking. Both builds treat compiler
warnings as errors and run all numerical tests and applications.

An `fpm.toml` manifest is supplied. `fpm` was not installed in the validation
environment, so the manifest is not claimed as tested.

## Applications

Fit a portfolio from a CSV file whose first column may be a date:

```sh
build/debug/bin/fit_csv data/example_returns.csv minvariance
build/debug/bin/fit_csv data/example_returns.csv tangency 0.0
build/debug/bin/fit_csv data/example_returns.csv riskparity
build/debug/bin/fit_csv data/example_returns.csv maxdiv
build/debug/bin/fit_csv data/example_returns.csv mad
build/debug/bin/fit_csv data/example_returns.csv cvar 0.05
build/debug/bin/fit_csv data/example_returns.csv efficient 0.0005
build/debug/bin/fit_csv data/example_returns.csv maxreturn
```

Run a rolling backtest:

```sh
build/debug/bin/backtest_csv data/example_returns.csv minvariance 60 20
```

`fit_csv` writes `weights.csv` in the current directory.

## Important numerical differences

The R package delegates optimization to numerous external backends. This
translation uses self-contained projected convex methods and active constraint
projections. The mean-variance solutions are tested against analytical values.
MAD uses smooth projected convex optimization. CVaR uses a primal-dual
solver for the exact finite-scenario Rockafellar-Uryasev dual formulation.
External GLPK solver paths and certificates are not claimed.

The tangency, diversification, and risk-parity routines use projected iterative
algorithms. Exact R optimizer endpoints are not claimed.

See `API_MAP.md`, `VALIDATION.md`, and `ORIGIN.md` for precise scope.
