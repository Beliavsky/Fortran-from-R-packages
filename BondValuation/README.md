# BondValuation for modern Fortran

A modern Fortran/FPM translation of the computational algorithms in the R package
`BondValuation` 0.1.1 by Wadim Djatschenko.

The library values fixed-coupon and zero-coupon bonds with regular or irregular
first and last coupon periods. It includes the package's complete calendar,
schedule, day-count, accrued-interest, price, yield, duration, and convexity
calculations.

## Main features

- Regular, short, and long first/last coupon periods
- Separate issue and first-interest-accrual dates
- Annual, semiannual, quarterly, triannual, bimonthly, and monthly coupons
- Zero-coupon bonds
- Sixteen day-count conventions
- Exact embedded Brazilian non-business-day calendar used by the R package
- Clean/dirty price conversion and accrued interest
- Clean price from yield and yield from clean price
- Modified and Macaulay duration
- Convexity using the original package convention
- Typed dates, schedules, bond terms, results, warning flags, and status codes
- Original-name compatibility module for all exported and low-level kernels

## Build with FPM

```sh
fpm build
fpm test
fpm run
fpm run --example regular_bond
fpm run --example day_count_comparison
```

The project has no external Fortran dependencies.

## Minimal example

```fortran
program example
  use bondvaluation
  implicit none
  type(bond_terms) :: terms
  type(bond_value_result) :: value

  terms%issue_date = date_from_ymd(2020, 1, 15)
  terms%maturity_date = date_from_ymd(2030, 1, 15)
  terms%coupon_frequency = 2
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.0_dp
  terms%day_count_convention = dcc_act_act_icma

  value = bond_price(4.0_dp, date_from_ymd(2024, 4, 15), terms)
  print *, value%clean_price
end program example
```

## Modules

- `bondvaluation`: preferred umbrella API
- `bondvaluation_dates`: Gregorian dates and schedule-date utilities
- `bondvaluation_daycount`: sixteen day-count conventions
- `bondvaluation_schedule`: coupon and anniversary schedules
- `bondvaluation_pricing`: accrued interest, price, yield, duration, convexity
- `bondvaluation_compat`: original-name and low-level compatibility API

Because Fortran is case-insensitive, the original compatibility procedure `DP`
would conflict with the common kind name `dp`. The compatibility module therefore
renames the real kind internally to `real_kind`; applications should similarly
rename one of the two when using both modules.

## Scope

R data frames, `timeDate` classes, Rcpp registration, formatted R warnings, and
bundled example-data presentation are not part of the numerical library. The
complete original package tree is retained in `original_source/` for provenance.
See `PORTING.md` for detailed behavior and differences.

## License

GPL-3.0-only, matching the original package metadata. See `LICENSE` and `NOTICE.md`.
