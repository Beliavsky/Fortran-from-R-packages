# Testing

The release test suite contains four independent programs.

## `test_series`

Checks return calculations and compounding, return conversion, capping, simple
and weighted moving averages, exponential smoothing, HP filtering, binary
trend conversion, matrix overloads, and matrix square roots.

## `test_risk`

Checks diversification and concentration ratios, volatility-weighted
correlation, the marginal-risk decomposition, and empirical/EVT tail dependence
for concordant and oppositely ranked series.

## `test_portfolios`

Checks global minimum-variance and most-diversified weights against independent
SciPy constrained-optimization references, verifies PMTD feasibility, and tests
equal-risk-contribution equality directly.

Reference weights used in the tests are:

- PGMV: `(0.14621071595666, 0.62650118429361, 0.22728809974973)`
- PMD: `(0.24110697, 0.45974267, 0.29915036)`

## `test_drawdown`

Checks all four drawdown formulations, hard and soft budgets, nonnegative
weights, drawdown feasibility, terminal-return references, and the zero-CDaR
solution of the deterministic test case. The terminal-return references were
also solved independently with `scipy.optimize.linprog`.

## Build configurations

Run:

```text
./run_tests.sh strict
./run_tests.sh optimized
```

The strict configuration enables warnings as errors, bounds and consistency
checks, floating-point traps, and backtraces. The optimized configuration uses
`-O3` and warnings as errors.
