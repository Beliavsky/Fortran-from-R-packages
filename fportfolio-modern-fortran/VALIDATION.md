# Validation

## Environment

- GNU Fortran 14.2.0
- LAPACK and BLAS from the Debian runtime environment
- Fortran 2018 source mode

## Debug flags

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace -ffree-line-length-none
```

## Release flags

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror
-fbacktrace -ffree-line-length-none
```

## Commands

```sh
./scripts/build_and_test.sh debug
./scripts/build_and_test.sh release
```

## Tests

`test_statistics_risk` checks:

- Sample and EWMA covariance symmetry
- Shrinkage covariance
- Covariance-risk Euler decomposition
- Historical, Normal, and modified VaR/ES relationships
- Historical and modified risk contributions and budgets
- Drawdown and maximum-loss calculations
- Empirical and Normal-margin tail-dependence matrices

`test_optimization` checks:

- Minimum-variance weights against an analytical diagonal-covariance solution
- Target-return equality
- Maximum-return corner solution
- Positive tangency Sharpe ratio
- Equal covariance-risk budgets
- Analytical diagonal maximum-diversification weights
- Efficient-frontier feasibility
- MAD and CVaR feasibility
- Exact cardinality constraints
- Generic linear inequality constraints

`test_backtest_monitor` checks:

- Rolling portfolio weights and budget constraints
- Transaction costs, turnover, wealth, and drawdown
- Rolling sigma, VaR, and CVaR
- MACD and drawdown indicators
- Turning points and rolling stability
- Net-performance deductions and rebalancing statistics

The build harness also runs every `fit_csv` optimizer mode and multiple
`backtest_csv` strategies.

## Scope of equivalence

The tests establish internal formulas, analytical special cases, feasibility,
and end-to-end execution. They do not establish iteration-by-iteration
identity with R's external optimizers. In particular, MAD is solved by a self-contained smooth projected method. CVaR is solved
by a self-contained primal-dual method rather than GLPK.
