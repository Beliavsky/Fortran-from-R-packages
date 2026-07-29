# API reference

Use the umbrella module for most programs:

```fortran
use nmof
```

All real-valued APIs use `real(dp)`. Procedures that can fail generally return an optional integer `status` or store it in a result object.

## Status codes

| Constant | Meaning |
|---|---|
| `nmof_ok` | Successful completion |
| `nmof_invalid_input` | Invalid dimensions, parameters, or option selector |
| `nmof_max_iter` | Iteration limit reached |
| `nmof_linear_solve_failed` | LAPACK or linear-system failure |
| `nmof_infeasible` | Constraints are infeasible |
| `nmof_not_bracketed` | Root was not bracketed |
| `nmof_numerical_failure` | Other numerical failure |

## Result types

| Type | Main fields |
|---|---|
| `optimization_result` | `xbest`, `ofvalue`, `population_values`, `history`, `iterations`, `status` |
| `binary_optimization_result` | logical `xbest` plus optimization diagnostics |
| `option_result` | `value`, standard Greeks, `status` |
| `cppi_result` | portfolio value, cushion, bond, floor, exposure, units, spot, `status` |
| `drawdown_summary` | maximum drawdown and high/low values and positions |
| `quadrature_rule` | `nodes`, `weights`, `status` |
| `frontier_result` | target returns, volatility, portfolio weights, `status` |
| `pbo_result` | PBO, logits, in-sample and out-of-sample scores, `status` |
| `qtable_result` | whiskers, median, minimum, maximum, `status` |
| `bond_return_result` | returns, duration, convexity, `status` |

## Optimization: `nmof_optimization`

| Procedure | Purpose |
|---|---|
| `de_opt` | Differential Evolution over bounded real vectors |
| `ps_opt` | Particle Swarm Optimization over bounded real vectors |
| `ga_opt` | Genetic Algorithm over logical/binary vectors |
| `local_search` | Iterative neighbourhood local search |
| `simulated_annealing` | Temperature-scheduled simulated annealing |
| `threshold_accepting` | Threshold Accepting with supplied or sampled thresholds |
| `greedy_search` | Best-improving search over an explicitly generated neighbourhood |
| `grid_search` | Cartesian product search over per-coordinate levels |
| `restart_opt` | Repeat a user optimizer and retain the best result |

Public callback interfaces:

- `real_objective`
- `real_neighbour`
- `real_repair`
- `real_velocity_change`
- `binary_objective`
- `binary_repair`
- `all_neighbours`
- `optimizer_callback`

Callbacks accept optional unlimited-polymorphic context objects, allowing typed user data without global variables.

## Portfolio methods: `nmof_portfolio`

| Procedure | Purpose |
|---|---|
| `minimum_variance` | Convex minimum-variance allocation |
| `mean_variance_portfolio` | Target-return or penalized mean-variance allocation |
| `mean_variance_frontier` | Multiple efficient-frontier portfolios |
| `maximum_sharpe` | Maximum excess-return-to-volatility portfolio |
| `tracking_portfolio` | Minimum tracking-error allocation |
| `equal_risk_contribution` | Equalize component variance contributions |
| `minimum_cvar` | Scenario-based CVaR minimization |
| `minimum_mad` | Scenario-based mean-absolute-deviation minimization |

The main convex portfolio routines support combinations of:

- full-investment budget constraints;
- per-asset lower and upper bounds;
- group membership matrices;
- group lower and upper bounds;
- minimum expected return.

## Fixed income and derivatives: `nmof_finance`

### Fixed income

| Procedure | Original NMOF counterpart |
|---|---|
| `vanilla_bond` | `vanillaBond` |
| `yield_to_maturity` | `ytm` with a scalar offset |
| `yield_to_maturity_curve` | `ytm` with vector offsets |
| `bond_duration` | `duration` |
| `bond_convexity` | `convexity` |
| `approximate_bond_return` | `approxBondReturn` |
| `bund_future` | `bundFuture` |
| `bund_future_implied_rate` | `bundFutureImpliedRate` |
| `xt_contract_value` | `xtContractValue` |
| `xt_tick_value` | `xtTickValue` |

### Options

| Procedure | Purpose |
|---|---|
| `vanilla_option_european` | European BSM call/put with optional Greeks and dividends |
| `vanilla_option_american` | American binomial call/put |
| `vanilla_option_implied_vol` | Recover volatility from European or American price |
| `put_call_parity` | Solve for a selected parity component |
| `barrier_option_european` | European up/down, in/out call/put barriers |
| `european_call_tree` | Direct recombining-tree European call |
| `european_call_binomial_expectation` | Closed binomial expectation form |

The variance argument is annualized variance, matching NMOF. Implied volatility is returned as volatility.

### Characteristic functions and Fourier pricing

| Procedure | Purpose |
|---|---|
| `cf_bsm` | Black-Scholes-Merton characteristic function |
| `cf_heston` | Heston stochastic-volatility characteristic function |
| `cf_bates` | Heston plus jumps |
| `cf_merton` | Merton jump-diffusion characteristic function |
| `cf_variance_gamma` | Variance-Gamma characteristic function |
| `call_cf` | Generic characteristic-function call pricer |
| `call_heston_cf` | Heston convenience call pricer |
| `call_merton` | Merton jump-diffusion call pricer |

`cf_callback` is the abstract interface for user-supplied characteristic functions.

## Simulation: `nmof_simulation`

| Procedure | Original counterpart |
|---|---|
| `geometric_brownian_motion` | `gbm` |
| `geometric_brownian_bridge` | `gbb` |
| `random_returns` | `randomReturns` |
| `resample_correlated` | `resampleC` |
| `cppi` | `CPPI` |

`random_returns` can enforce requested finite-sample means, standard deviations, and correlation exactly, subject to feasible dimensions and a positive-semidefinite target matrix.

## Curves, risk, and general utilities: `nmof_utilities`

| Procedure | Purpose |
|---|---|
| `moving_average` | NMOF moving average |
| `partial_moment` | Lower or upper partial moment |
| `conditional_moment` | Conditional lower or upper moment |
| `drawdown_series` | Absolute or relative drawdown path |
| `drawdown_info` | Maximum drawdown summary |
| `diversification_ratio` | Weighted asset volatility divided by portfolio volatility |
| `ns_curve` | Nelson-Siegel curve |
| `nss_curve` | Nelson-Siegel-Svensson curve |
| `ns_factors` | Nelson-Siegel factor matrix |
| `nss_factors` | Nelson-Siegel-Svensson factor matrix |
| `change_interval` | Transform quadrature nodes and weights to another interval |
| `xw_gauss` | Gaussian quadrature nodes and weights |
| `bracketing` | Locate sign-change intervals |
| `repair_matrix` | Eigenvalue-clipped symmetric matrix repair |
| `column_subset` | Rank-check candidate column subsets |
| `qtable_statistics` | Quantile-table summary |
| `probability_backtest_overfitting` | Combinatorially symmetric cross-validation PBO |
| `marginal_risk_contributions` | Analytical covariance-based MRC |
| `marginal_risk_contributions_fd` | General finite-difference MRC callback |

Public risk callback interfaces:

- `metric_function`
- `scalar_function_context`

### Test functions

- `test_ackley`
- `test_eggholder`
- `test_griewank`
- `test_rastrigin`
- `test_rosenbrock`
- `test_schwefel`
- `test_trefethen`

## Numerical support APIs

### `nmof_rng`

- `rng_state`
- `rng_seed`
- `rng_uniform`
- `rng_normal`
- `rng_integer`
- `rng_logical`
- `rng_shuffle`

### `nmof_math`

- normal PDF, CDF, and quantile;
- R type-7 quantile, median, standard deviation, and sorting;
- bisection and Brent root solving;
- adaptive Simpson integration;
- factorial/log-factorial helpers and clamping.

### `nmof_linalg`

- general and symmetric-positive-definite solves;
- matrix inversion;
- eigendecomposition and Cholesky factorization;
- covariance matrices, column means/standard deviations, and rank tests.

### `nmof_qp`

- equality-constrained quadratic programming;
- active-set quadratic programming;
- equality, feasible-set, linear-constraint, and budget-box projection;
- feasibility checks.

## Original exported-name mapping

| R export | Fortran status |
|---|---|
| `CPPI` | `cppi` |
| `DEopt` | `de_opt` |
| `EuropeanCall` | `european_call_tree` |
| `EuropeanCallBE` | `european_call_binomial_expectation` |
| `GAopt` | `ga_opt` |
| `LSopt` | `local_search` |
| `MA` | `moving_average` |
| `NS`, `NSS` | `ns_curve`, `nss_curve` |
| `NSf`, `NSSf` | `ns_factors`, `nss_factors` |
| `PBO` | `probability_backtest_overfitting` |
| `PSopt` | `ps_opt` |
| `SAopt` | `simulated_annealing` |
| `TAopt` | `threshold_accepting` |
| `approxBondReturn` | `approximate_bond_return` |
| `barrierOptionEuropean` | `barrier_option_european` |
| `bracketing` | `bracketing` |
| `bundFuture` | `bund_future` |
| `bundFutureImpliedRate` | `bund_future_implied_rate` |
| `callCF` | `call_cf` |
| `callHestoncf` | `call_heston_cf` |
| `callMerton` | `call_merton` |
| `cfBSM`, `cfBates`, `cfHeston`, `cfMerton`, `cfVG` | corresponding `cf_*` procedures |
| `changeInterval` | `change_interval` |
| `colSubset` | `column_subset` |
| `convexity` | `bond_convexity` |
| `divRatio` | `diversification_ratio` |
| `drawdown` | `drawdown_series`, `drawdown_info` |
| `duration` | `bond_duration` |
| `gbb`, `gbm` | geometric Brownian bridge/motion |
| `greedySearch`, `gridSearch` | `greedy_search`, `grid_search` |
| `maxSharpe` | `maximum_sharpe` |
| `minCVaR`, `minMAD`, `minvar` | corresponding portfolio procedures |
| `mvFrontier`, `mvPortfolio` | corresponding mean-variance procedures |
| `pm` | `partial_moment`, `conditional_moment` |
| `putCallParity` | `put_call_parity` |
| `qTable` | `qtable_statistics` |
| `randomReturns` | `random_returns` |
| `repairMatrix` | `repair_matrix` |
| `resampleC` | `resample_correlated` |
| `restartOpt` | `restart_opt` |
| `tf*` | `test_*` benchmark functions |
| `trackingPortfolio` | `tracking_portfolio` |
| `vanillaBond` | `vanilla_bond` |
| `vanillaOption*` | corresponding `vanilla_option_*` procedures |
| `xtContractValue`, `xtTickValue` | corresponding XT procedures |
| `xwGauss` | `xw_gauss` |
| `ytm` | `yield_to_maturity`, `yield_to_maturity_curve` |
| `French`, `Ritter`, `Shiller` | omitted network/data acquisition |
| `showChapterNames`, `showExample` | omitted R package discovery |
| `LS.info`, `SA.info`, `TA.info` | replaced by callback iteration arguments/result history |
