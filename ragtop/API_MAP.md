# API map

## Direct numerical counterparts

| R API | Fortran API |
|---|---|
| `blackscholes` | `black_scholes` |
| `black_scholes_on_term_structures` | `black_scholes_on_term_structures` |
| `variance_cumulation_from_vols` | `initialize_volatility_curve`, `cumulative_variance` |
| `spot_to_df_fcn` | `initialize_discount_curve`, `discount_factor` |
| `construct_implicit_grid_structure` | same name |
| `construct_tridiagonals` | same name |
| `pde_matrix_solve` | same name |
| `infer_conforming_time_grid` | same name |
| `form_present_value_grid` | same name |
| `find_present_value` | same name |
| `american` | same name |
| `implied_volatility` | same name |
| `implied_volatilities` | same name |
| `implied_volatility_with_term_struct` | same name |
| `american_implied_volatility` | same name |
| `equivalent_jump_vola_to_bs` | same name |
| `equivalent_bs_vola_to_jump` | same name |
| `implied_jump_process_volatility` | same name |
| `fit_variance_cumulation` | same name |
| `price_with_intensity_link` | same name |
| `penalty_with_intensity_link` | same name |
| `fit_to_option_market` | same name |
| `shift_for_dividends` | same name |
| `time_adj_dividends` | `time_adjusted_dividends` |
| `adjust_for_dividends` | same name |
| `value_from_prior_coupons` | same name |
| `accelerated_coupon_value` | same name |
| `coupon_value_at_exercise` | same name |
| `construct_descending_bumps` | same name |
| `resolve_bumps` | `resolve_bump` |
| `greek_by_fd` | same name |
| `robust_greek` | same name |
| `grid_delta_gamma` | same name |
| `find_greeks` | same name |
| `check_discount_factor_fcn` | `check_discount_factor` |
| `check_variance_cumulation_fcn` | `check_variance_cumulation` |
| `check_survival_probability_fcn` | `check_survival_probability` |

## Instrument classes

The R reference classes are represented by a single typed `instrument_spec`
with constructor functions:

- `EuropeanOption`
- `AmericanOption`
- `ZeroCouponBond`
- `CouponBond`
- `CallableBond`
- `ConvertibleBond`

`terminal_values`, `recovery_values`, `apply_optionality`, and
`instrument_cashflow_between` implement the computational class behavior.

## Adapted interfaces

- R function closures for rates, volatility, and hazard are represented by
  typed `market_spec`, `discount_curve`, `volatility_curve`, and `hazard_spec`.
- `fit_to_option_market` uses deterministic bounded coordinate search instead
  of optional R optimization packages.
- Callable and putable schedules are active in the Fortran implementation;
  the upstream `CallableBond` class records those schedules but does not
  override its optionality method in the supplied source.
- Multiple instruments are priced by calling `find_present_value` for each
  typed instrument. The original shared R grid was primarily an efficiency and
  state-sharing mechanism.

## Omitted R-specific interfaces

- Reference-class and S4 metadata, `$` methods, caches, and logging
- Data-frame conversion wrapper `fit_to_option_market_df`
- `treasury_df_raw` and `treasury_df`, which download or cache external data
- `detail_from_AnnivDates`, which depends on R date objects and
  `BondValuation`
- Plotting, vignettes, packaged `.rda` data loading, and optional R packages
