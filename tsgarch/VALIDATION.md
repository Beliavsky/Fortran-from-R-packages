# Validation

Six deterministic test programs cover the translated computational surface.

## `test_models`

- Filters all eight volatility-model families.
- Verifies strictly positive conditional variances.
- Evaluates all ten innovation distributions inside a GARCH filter.
- Confirms that equation descriptions are available.

## `test_simulation_forecast`

- Runs unconditional multi-path Student-t GARCH simulation.
- Runs history-conditioned simulation.
- Checks dimensions, success statuses, and positive forecast sigmas.
- Verifies ordered 5%, 50%, and 95% forecast quantiles.

## `test_estimation`

- Simulates a stationary Gaussian GARCH process.
- Fits the model from generic initialized parameters.
- Verifies that optimization improves the log likelihood.
- Checks nonnegative coefficients and fitted stationarity.

## `test_inference_profile`

- Constructs numerical Hessian and per-observation score matrices.
- Checks PIT residuals lie in the unit interval.
- Computes Wald confidence intervals and OPG covariance.
- Evaluates a likelihood profile and persistence half-life.

## `test_backtest`

- Runs a rolling one-step 5% VaR backtest.
- Checks forecast dimensions and refit count.
- Checks empirical coverage and Kupiec p-values are valid probabilities.
- The returned object also contains Christoffersen independence and conditional-
  coverage statistics.

## `test_constraints_vreg`

- Filters a variance-targeted model with an external variance regressor.
- Verifies parameter pack/unpack round trips with regressor coefficients.
- Verifies the IGARCH coefficient equality to machine precision.

## Build modes

The checked GNU Fortran build uses:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic
-fcheck=all -fbacktrace -O0
```

The optimized build uses:

```text
-std=f2018 -Wall -Wextra -Werror -pedantic -O3 -march=native
```

Both modes compile and pass all six tests without warnings.

## Demonstration

The example simulates a Student-t GJR-GARCH process, estimates it, and produces
five-step Monte Carlo volatility forecasts. Its purpose is an end-to-end API
smoke test rather than a fixed statistical benchmark.

## Expected cross-language differences

- Fortran and R use different random-number streams.
- The bounded Nelder-Mead optimizer differs from TMB/nloptr.
- Hessians and scores use finite differences instead of automatic
  differentiation.
- Forecasts use Monte Carlo for all models.

Consequently, stochastic output and estimated final digits are not expected to
match R bit for bit. Structural constraints, likelihood improvement, recursion
identities, dimensions, and probability ranges are tested directly.
