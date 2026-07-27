# rtl-fortran

Modern Fortran 2018 and FPM translation of the numerical algorithms in RTL 1.3.9.

The original R package is a commodities risk and analytics toolkit. This port focuses on its self-contained numerical routines and intentionally excludes GARCH, as requested. It also excludes web APIs, plotting, R object infrastructure, and bundled market-data loading.

## License

MIT, matching the supplied RTL source package. See `LICENSE`, `NOTICE`, and the retained upstream source under `original/RTL-1.3.9`.

## Build

```text
fpm build
fpm test
fpm run rtl_demo
fpm run --example options_and_processes
fpm run --example portfolio_and_trading
fpm run --example upstream_names
```

The release contains 2,703 lines across 19 Fortran source, test, demo, and example files.

The library has no external dependencies.

## Main modules

- `rtl_options`: generalized Black-Scholes, CRR trees, Kirk spread options, and barrier spread options.
- `rtl_processes`: GBM, OU, time-varying OU, OU with jumps, OU fitting, and multivariate simulation.
- `rtl_fixed_income`: bonds, NPV, interest-rate swaps, and commodity swaps.
- `rtl_calendar`: Gregorian dates, business-day counting, and commodity futures weights.
- `rtl_portfolio`: random portfolio clouds, a native simplex solver, and refinery optimization.
- `rtl_market`: returns, roll masks, prompt betas, performance statistics, and moving-average strategies.
- `rtl`: umbrella module exporting the public API.

## Example

```fortran
program example
  use rtl, only: dp, option_result, gbs_option
  implicit none

  type(option_result) :: result

  result = gbs_option(100.0_dp, 100.0_dp, 1.0_dp, 0.05_dp, &
    0.02_dp, 0.20_dp, "call")

  print '(f12.6)', result%price
end program example
```

Compatibility-style names such as `GBSOption`, `CRROption`, `simOU`, `fitOU`, `swapIRS`, `promptBeta`, and `tradeStrategySMA` are provided in addition to descriptive snake-case names.

## Deliberate design choices

- `sim_multivariates` honors a supplied starting vector. Set `use_last_start=.true.` to reproduce the upstream behavior that replaces it with the final historical observation.
- `sim_ou_jump` draws independent jump sizes by default. Set `legacy_recycled_jump=.true.` to reproduce the upstream single recycled jump size.
- `npv_value` preserves the upstream convention that terminal value replaces the final periodic cash flow. Set `terminal_replaces_cash_flow=.false.` to add it instead.
- `interest_rate_swap` supports semiannual resets, which are documented upstream but omitted from its branch logic.
- `barrier_spread_option` accepts `monitoring` for compatibility, but reports `monitoring_used=.false.` because the upstream formula does not use it.
- `efficient_frontier` preserves the upstream random portfolio-cloud method; it is not presented as an exact constrained frontier optimizer.
- `trade_stats` supplies dependency-free equivalents of the upstream external package metrics. See `PORTING_NOTES.md` for formula details.

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for the complete mapping and validation record.
