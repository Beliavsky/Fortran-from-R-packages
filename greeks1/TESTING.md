# Testing

## FPM

```text
fpm test
```

## Direct GNU Fortran

On Linux or macOS:

```text
./scripts/test_gnu.sh
```

On Windows with GNU Fortran on `PATH`:

```text
scripts\test_gnu.bat
```

The direct build uses:

```text
-std=f2018 -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

## Test programs

- `test_black_scholes`: known Black-Scholes values, Greeks, parity, and binary
  options.
- `test_geometric_binomial`: upstream geometric-Asian reference values and
  American-option tree behavior.
- `test_implied_helpers`: implied-volatility recovery, dispatcher behavior,
  Brownian path layout, and quadrature helpers.
- `test_monte_carlo`: European, geometric-Asian, arithmetic-Asian,
  control-variate, and jump-diffusion Monte Carlo paths.

Monte Carlo tests use fixed seeds and statistical tolerances rather than exact
path matching with R's different random-number engine.
