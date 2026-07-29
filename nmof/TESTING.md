# Testing

## Test commands

```sh
./scripts/run_tests.sh strict
./scripts/run_tests.sh optimized
```

The strict configuration uses GNU Fortran with:

```text
-O0 -g -std=f2018 -ffree-line-length-none
-Wall -Wextra -Werror
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The optimized configuration uses `-O3` and rebuilds every source, test, application, and example from a clean directory.

## Test programs

### `test_core`

Covers:

- moving averages and partial moments;
- drawdown calculations;
- benchmark objective functions;
- Gauss quadrature and transformed intervals;
- root bracketing;
- positive-semidefinite matrix repair;
- minimum-variance and ERC identities;
- exact sample means, standard deviations, and correlations from `random_returns`;
- analytical and finite-difference marginal risk contributions;
- correlated resampling;
- GBM and Brownian bridges;
- CPPI accounting identities;
- Probability of Backtest Overfitting.

### `test_finance`

Covers:

- vanilla bond pricing and yield recovery;
- duration, convexity, and approximate returns;
- European Black-Scholes-Merton prices and Greeks;
- discrete dividends;
- implied-volatility recovery;
- American options;
- binomial European-call identities;
- put-call parity;
- barrier options;
- characteristic functions and Heston Fourier pricing.

### `test_optimization`

Exercises:

- Differential Evolution;
- Particle Swarm Optimization;
- binary Genetic Algorithms;
- Local Search;
- Simulated Annealing;
- Threshold Accepting;
- Greedy Search;
- Grid Search;
- Restart Optimization.

The tests use deterministic seeds and objectives with known optima or independently checked results.

### `test_portfolio`

Covers:

- bounded and group-constrained minimum variance;
- mean-variance target-return constraints;
- efficient frontiers;
- maximum Sharpe ratio;
- tracking-portfolio recovery;
- minimum-CVaR and minimum-MAD budget/bound constraints.

## Reference sources

Reference values come from a combination of:

- formulas and examples in the original NMOF source and documentation;
- the original `inst/tinytest` and `inst/unitTests` trees retained under `original/`;
- independent analytical identities;
- independently assembled covariance, optimization, and quadrature calculations.

The original R package was not required at runtime for the Fortran build.

## Linker note

GNU ld may report that objects containing internal callback procedures require an executable stack. This is caused by GNU Fortran trampolines for nested procedures used as callbacks; it is a linker warning rather than a numerical test failure.
