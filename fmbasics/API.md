# API reference

All public entities are available through `use fmbasics`.

## Status and precision

- `dp`: double-precision real kind.
- `FM_OK`: successful operation.
- `FM_INVALID_ARGUMENT`, `FM_DOMAIN_ERROR`, `FM_IO_ERROR`,
  `FM_NOT_CONVERGED`: nonzero status values.

Most constructors and numerical routines accept an optional final `status`
argument. Functions return an empty or zero-valued object when validation
fails; callers requiring diagnostics should inspect `status`.

## Dates and calendars

- `make_date(year, month, day)`
- `date_from_yyyymmdd(value)`, `date_to_yyyymmdd(date)`
- `date_year`, `date_month`, `date_day`, `date_string`
- `days_period`, `months_period`, `years_period`
- `calendar`, `joint_calendar`, `calendar_contains`
- `is_good_day`, `shift_date`, `roll_date`, `add_months`
- `year_frac(from_date, to_date, basis)`

Supported day-count names include `act/365`, `act/360`, `act/act`,
`30e/360` and `30/360us`.

## Currencies, FX and indices

Types:

- `currency_t`
- `currency_pair_t`
- `index_t`

Currency constructors: `aud`, `eur`, `gbp`, `jpy`, `nzd`, `usd`, `chf`,
`hkd`, `nok`.

Currency-pair constructors: `audusd`, `eurusd`, `nzdusd`, `gbpusd`,
`usdjpy`, `gbpjpy`, `eurgbp`, `audnzd`, `eurchf`, `usdchf`, `usdhkd`,
`eurnok`, `usdnok`.

Operations:

- `pair_iso(pair)`
- `invert(pair)`
- `is_t1(pair)`
- `to_spot`, `to_spot_next`, `to_today`, `to_tomorrow`
- `to_forward(date, tenor, pair)`
- `to_fx_value(date, name_or_period, pair)`

Benchmark term indices:

- `audbbsw`, `audbbsw1b`, `euribor`, `gbplibor`, `jpylibor`, `jpytibor`
- `nzdbkbm`, `usdlibor`, `chflibor`, `hkdhibor`, `noknibor`

Benchmark cash indices:

- `aonia`, `eonia`, `sonia`, `tonar`, `nziona`, `fedfunds`, `chftois`,
  `honix`

Index-date functions: `to_reset`, `to_value`, `to_maturity`.

## Rates and discount factors

Types:

- `interest_rate_t`
- `discount_factor_t`

Constructors and conversions:

- `interest_rate(value, compounding, day_basis)`
- `discount_factor(value, start_date, end_date)`
- `as_discount_factor(rate, start_date, end_date)`
- `as_interest_rate(df, compounding, day_basis)`
- `convert_interest_rate(rate, compounding, day_basis)`
- `compound_factor`, `implied_rate`, `is_valid_compounding`

`COMPOUND_CONTINUOUS` denotes continuous compounding. Numeric compounding
values such as 1, 2, 4, 12 and 365 denote the number of compounds per year;
0 denotes simple interest and -1 denotes simple discount.

Explicit arithmetic functions are `rate_add`, `rate_subtract`,
`rate_multiply`, `rate_divide`, `discount_multiply` and `discount_divide`.

## Interpolation and zero curves

Types:

- `interpolation_t`
- `zero_curve_t`

Interpolation constructors:

- `constant_interpolation`
- `linear_interpolation`
- `cubic_interpolation`
- `logdf_interpolation`
- `linear_cubic_time_var_interpolation`

Zero-curve functions:

- `zero_curve(discount_factors, reference_date, interpolation)`
- `build_zero_curve(interpolation, filename)`
- `load_zero_curve_csv`
- `interpolate_zero(curve, terms)`
- `interpolate_zeros(curve, dates)`
- `interpolate_dfs(curve, from_dates, to_dates)`
- `interpolate_fwds(curve, from_dates, to_dates)`

## Credit

Types:

- `cds_spec_t`, `cds_curve_t`
- `survival_probabilities_t`, `zero_hazard_rate_t`
- `credit_curve_t`

Constructors and conversions:

- `cds_spec`, `cds_single_name_spec`, `cds_markit_spec`, `cds_curve`
- `survival_probabilities`, `zero_hazard_rate`
- `as_survival_probabilities`, `as_zero_hazard_rate`
- `survival_multiply`, `survival_divide`
- `bootstrap_cds_survival`
- `credit_curve`
- `interpolate_credit`, `interpolate_credit_zeros`
- `interpolate_credit_dfs`, `interpolate_credit_fwds`

## Volatility surfaces

Types:

- `vol_quotes_t`
- `vol_surface_t`

Functions:

- `vol_quotes`, `vol_surface`
- `load_vol_quotes_csv`, `build_vol_quotes`, `build_vol_surface`
- `interpolate_vol(surface, maturity_dates, smile_coordinates)`

Interpolation is linear in total variance across maturity and natural cubic
across the smile coordinate, matching the original package method.

## Money and cash flows

Types:

- `single_currency_money_t`
- `multi_currency_money_t`
- `cash_flow_t`

Functions:

- `single_currency_money`
- `multi_currency_money`
- `combine_money`
- `aggregate_by_currency`
- `cash_flow`

The derived types expose `%size()` methods. Values and currencies are stored in
parallel allocatable arrays for straightforward numerical use.
