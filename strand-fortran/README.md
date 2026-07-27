# strand-fortran

A modern Fortran 2018/FPM numerical port of the portfolio-construction and
share-level simulation engine in [`strand` 0.2.3](original/strand-0.2.3).

The upstream package is an R6-based framework for discrete investment-strategy
simulation. This port extracts its array-oriented numerical algorithms into a
self-contained Fortran library. It preserves the upstream GPL version 3 license
and attribution.

## Implemented numerical scope

- rank-normal signal transformation, grouped normalization, and repeated
  numeric-factor neutralization;
- cross-section carry-forward, missing-value replacement, and period statistics;
- internal/external share bookkeeping and corporate-action adjustment ratios;
- multi-strategy linear-program portfolio construction;
- long/short target market values, position and trading limits, turnover limits,
  factor constraints, category constraints, and constraint loosening;
- exact joint-order netting into internal transfers and external market orders;
- volume-limited fills, transaction costs, financing costs, dividends,
  distributions, delistings, and daily P&L;
- multi-day share-level simulation;
- factor/category exposure aggregation and performance statistics;
- a dependency-free two-phase simplex LP solver.

See [COVERAGE.md](COVERAGE.md) for the detailed mapping and
[PORTING_NOTES.md](PORTING_NOTES.md) for behavioral and API differences.

## Build with FPM

```text
fpm build
fpm test
fpm run strand_demo
fpm run --example portfolio_optimization
fpm run --example share_level_simulation
```

The project has no external numerical-library dependency.

## Minimal portfolio optimization

```fortran
program example
  use strand
  implicit none
  type(optimizer_config) :: config
  type(optimization_result) :: result
  real(dp) :: price(4), dollar_volume(4), alpha(4, 1)
  logical :: investable(4)
  integer :: shares(4, 1)

  allocate(config%strategies(1))
  config%strategies(1)%capital = 1000000.0_dp
  config%strategies(1)%ideal_long_weight = 0.5_dp
  config%strategies(1)%ideal_short_weight = 0.5_dp
  config%strategies(1)%target_long_weight = 0.5_dp
  config%strategies(1)%target_short_weight = 0.5_dp
  config%strategies(1)%has_target_weights = .true.
  config%strategies(1)%position_limit_pct_adv = 100.0_dp
  config%strategies(1)%position_limit_pct_lmv = 100.0_dp
  config%strategies(1)%position_limit_pct_smv = 100.0_dp
  config%strategies(1)%trading_limit_pct_adv = 100.0_dp

  price = 50.0_dp
  dollar_volume = 10000000.0_dp
  alpha(:, 1) = [2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp]
  investable = .true.
  shares = 0

  result = optimize_portfolio(config, price, dollar_volume, investable, alpha, shares)
  if (.not. result%success) error stop trim(result%message)
  print *, result%order_shares(:, 1)
end program example
```

## Main modules

| Module | Purpose |
|---|---|
| `strand` | Public umbrella module |
| `strand_types` | Typed configurations and results |
| `strand_stats` | normalization, neutralization, exposures, performance |
| `strand_data` | cross-section and portfolio state |
| `strand_optimizer` | constrained multi-strategy portfolio construction |
| `strand_simulation` | fills, transfers, costs, P&L, and lifecycle simulation |
| `strand_simplex` | native bounded LP solver |
| `strand_linalg` | least squares and correlation support |

## R-specific functionality not compiled

YAML parsing, R6 classes, data frames, Arrow/Feather storage, Shiny, plotting,
report generation, tidyverse adapters, and bundled sample-data ingestion remain
in the retained original source but are not part of the Fortran library.
Callers provide already aligned numeric arrays and typed configuration objects.

## Validation

Run the direct compiler validation when FPM is not installed:

```text
./scripts/validate.sh
```

On Windows:

```bat
scripts\validate.bat
```

The validation builds checked and optimized configurations, runs all tests,
and executes the demo and examples. See [VALIDATION.md](VALIDATION.md).

## License

`strand-fortran` is licensed under **GPL-3.0-only**, matching the license declared
by `strand` 0.2.3. See [LICENSE](LICENSE) and [NOTICE](NOTICE).
