# API map

## Direct numerical counterparts

| R function or class | Fortran counterpart |
|---|---|
| `calc_gompertz_survival_probability()` | `gompertz_survival_probability()` |
| internal Gompertz density | `gompertz_pdf()` |
| `calc_gompertz_parameters()` | `fit_gompertz_mortality()` |
| `calc_gompertz_joint_parameters()` | `fit_joint_gompertz()` and `household_joint_survival()` |
| `calc_life_expectancy()` | `life_expectancy()` |
| `calc_gompertz_mode()` | `gompertz_mode_from_life_expectancy()` |
| `calc_retirement_ruin()` | `retirement_ruin_probability()` |
| internal `calc_a()` | `gompertz_annuity_factor()` |
| internal incomplete gamma | `upper_incomplete_gamma()` |
| `calc_purchasing_power()` | `purchasing_power()` |
| `calc_optimal_risky_asset_allocation()` | `optimal_risky_asset_allocation()` |
| `calc_risk_adjusted_return()` | `risk_adjusted_return()` |
| internal `calc_present_value()` | `present_value_stream()` |
| internal utility/inverse utility | `utility()` and `inverse_utility()` |
| internal certainty-equivalent return | `certainty_equivalent_return()` |
| `calc_effective_tax_rate()` | `effective_tax_rates()` |
| `calc_portfolio_parameters()` | `calculate_portfolio_parameters()` |
| internal portfolio covariance/moments | `covariance_from_sd_corr()`, `calculate_joint_networth_moments()` |
| internal `calc_expected_utility()` | `expected_utility()` |
| internal net-worth fractions | `networth_fractions()` |
| `calc_optimal_portfolio()` | `optimize_portfolio()` |
| `create_portfolio_template()` | `create_default_portfolio()` |
| internal random-return generator | `generate_random_returns()` |
| `HouseholdMember` | `household_member` |
| `Household` | `household` |
| internal household timeline | `build_household_timeline()` |
| internal cash-flow generation | `generate_cashflow_stream()` with a procedure callback |
| internal discretionary spending | `discretionary_spending()` |
| `simulate_single_scenario()` | `simulate_lifecycle()` |
| Monte Carlo portion of `simulate_scenario()` | `simulate_lifecycle_samples()` |

## Adapted interfaces

- R expressions such as `"age <= 65 ~ 70000"` are not interpreted. A Fortran
  callback receives the period, date, and member ages and returns a cash flow.
- `simulate_lifecycle()` receives already assembled income and essential
  spending vectors. This separates financial calculations from user-interface
  rule parsing.
- `optimize_portfolio()` uses an analytic-gradient projected method on one or
  two simplices rather than `nloptr::NLOPT_LD_SLSQP`. It reproduces the
  upstream benchmark allocations to the shown precision.
- Joint household survival is fitted with native bounded Nelder-Mead.
- Calendar ages use Gregorian day counts divided by 365.2425.

## Not translated

- Shiny application and UI helpers
- all `plot_*` and `render_scenario_snapshot()` routines
- R6/S3 printing, nested tibbles, and tidy-evaluation infrastructure
- expression parsing in `generate_cashflow_streams()`
- caching, progress bars, `future`/`furrr` parallelism
- bundled `.rda` HMD life tables and `get_default_gompertz_parameters()`
- HMD text-file parsing and package data preparation
- pkgdown, spelling, snapshots, and visual regression tests

Users can pass externally loaded age/mortality vectors to
`fit_gompertz_mortality()` and can run independent simulation samples in
parallel at the application level.
