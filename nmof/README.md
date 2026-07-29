# NMOF for Modern Fortran

This project is an independent modern Fortran translation of the computational algorithms in **NMOF 2.12-0**, *Numerical Methods and Optimization in Finance*.

The code uses a standard Fortran Package Manager layout, Fortran 2018, explicit kinds, `implicit none`, typed result objects, callback interfaces, deterministic random-number streams, and explicit status codes. BLAS and LAPACK are the only external numerical libraries.

## Implemented computational areas

### Optimization

- Differential Evolution
- Particle Swarm Optimization
- Genetic Algorithms for binary decision vectors
- Local Search
- Simulated Annealing
- Threshold Accepting
- Greedy Search
- Grid Search
- Repeated/restarted optimization
- Callback-based objectives, neighbourhoods, repairs, penalties, velocity changes, and iteration reporting
- Ackley, Eggholder, Griewank, Rastrigin, Rosenbrock, Schwefel, and Trefethen test functions

### Portfolio methods and risk

- Minimum-variance portfolios
- Mean-variance portfolios and efficient frontiers
- Maximum-Sharpe portfolios
- Tracking portfolios
- Equal-risk-contribution portfolios
- Minimum-CVaR portfolios
- Minimum-MAD portfolios
- Asset bounds, budget constraints, group constraints, and minimum-return constraints
- Marginal risk contributions, diversification ratio, partial moments, and conditional moments
- Drawdown series and maximum-drawdown summaries
- Probability of Backtest Overfitting (PBO)

### Fixed income and term structures

- Vanilla bond pricing
- Yield-to-maturity solving
- Duration and convexity
- Approximate bond total returns
- Nelson-Siegel and Nelson-Siegel-Svensson curves and factor matrices
- German Bund futures and implied-rate calculations
- XT contract and tick values

### Options and characteristic-function pricing

- European Black-Scholes-Merton calls and puts, with Greeks
- American binomial calls and puts
- Implied volatility
- Put-call parity
- European barrier options
- Direct binomial European-call valuation
- Black-Scholes-Merton, Heston, Bates, Merton, and Variance-Gamma characteristic functions
- Fourier call pricing from arbitrary characteristic-function callbacks
- Heston and Merton convenience pricers

### Simulation and numerical utilities

- Geometric Brownian motion and Brownian bridges
- Exact-sample-moment random return generation
- Correlated resampling
- Constant Proportion Portfolio Insurance
- Moving averages
- Gaussian quadrature and interval transformation
- Root bracketing, bisection, Brent solving, and adaptive Simpson integration
- Covariance/correlation matrix repair
- Quantile-table summaries and column-subset rank checks
- Portable random-number generation

## Build with FPM

```sh
fpm build
fpm test
fpm run
fpm run --example finance_example
fpm run --example optimization_example
```

The package links to LAPACK and BLAS as declared in `fpm.toml`.

## Reproducible GNU Fortran build

When FPM is unavailable, the included script builds the modules in dependency order and runs all tests, applications, and examples:

```sh
./scripts/run_tests.sh strict
./scripts/run_tests.sh optimized
```

## Minimal example

```fortran
program example
   use nmof, only: dp, option_result, vanilla_option_european
   implicit none

   type(option_result) :: option

   option = vanilla_option_european( &
      spot=100.0_dp, strike=100.0_dp, tau=1.0_dp, &
      r=0.03_dp, q=0.01_dp, variance=0.20_dp**2, &
      option_type='call', greeks=.true.)

   print '(a,f12.6)', 'call value = ', option%value
   print '(a,f12.6)', 'delta      = ', option%delta
end program example
```

NMOF follows the original convention that the option routines receive **variance**, not volatility.

## Scope differences

Network downloads (`French`, `Ritter`, and `Shiller`), plotting, S3 printing, R call-stack inspection, progress bars, parallel orchestration, dataset loading, and example-file discovery are R runtime or presentation facilities rather than numerical algorithms and are not reproduced. The original package tree is retained under `original/` for provenance.

Detailed function mapping and intentional numerical differences are described in `API.md` and `PORTING.md`.

## License and provenance

The original NMOF package is Copyright (C) 2010-2025 Enrico Schumann and licensed under GPL-3. This translation preserves the **GPL-3.0-only** license. See `LICENSE`, `NOTICE.md`, and the retained `original/` source tree.

This is an independent translation and is not an official release of the original NMOF project.
