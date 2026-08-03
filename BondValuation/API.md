# API reference

## Preferred umbrella module

```fortran
use bondvaluation
```

### Core types

`date_type(year, month, day)` stores a Gregorian civil date.

`bond_terms` contains:

- `issue_date`, `maturity_date`
- `coupon_frequency`: `0`, `1`, `2`, `3`, `4`, `6`, or `12`
- optional `first_coupon_date`, `last_coupon_date`
- optional `first_interest_accrual_date`
- `redemption_value`, `coupon_rate_percent`
- `day_count_convention`
- `end_of_month`, `infer_end_of_month`
- `regular_cashflows_equal`

`bond_schedule` returns real cash-flow dates, regular anniversary dates, period
coordinates, coupon amounts, odd-period classifications, inferred dates, warnings,
and a status code.

`bond_value_result` returns clean and dirty prices, accrued interest, yield,
modified/Macaulay duration in years and coupon periods, convexity, settlement
coordinate, iteration count, and status.

### Dates

- `date_from_ymd(year, month, day)`
- `valid_date(date)`
- `leap_year(year)`
- `days_in_month(year, month)` and `days_in_year(year)`
- `date_to_serial(date)` and `serial_to_date(serial)` using R's 1970-01-01 origin
- `day_diff(first, second)`
- `add_months(date, months, end_of_month, reference_day)`
- `last_day_of_month(date)`, `is_last_day_of_month(date)`
- `date_to_string(date)`
- sorting and previous/next-date helpers

### Day-count conventions

Constants `dcc_act_act_isda` through `dcc_bus_252` correspond to:

1. Actual/Actual ISDA
2. Actual/Actual ICMA
3. Actual/Actual AFB
4. Actual/365L
5. 30/360 Bond Basis
6. 30E/360
7. 30E/360 ISDA
8. 30/360 German
9. 30U/360 US
10. Actual/365 Fixed
11. Actual (No Leap)/365
12. Actual/360
13. 30/365
14. Actual/365 Canadian Bond
15. Actual/364
16. Business/252 Brazil

```fortran
result = day_count_fraction(first_date, second_date, dcc, &
  coupon_frequency, maturity, end_of_month, next_coupon_year, anniversary_dates)
```

`result` is a `daycount_result` containing `days_accrued`, `fraction`, and `status`.
`year_fraction` returns only the fraction. `day_count_name` returns a label.

### Schedule generation

```fortran
call build_bond_schedule(terms, schedule)
call anniv_dates(terms, schedule)
```

Related routines are `period_coordinate`, `previous_coupon_date`, and
`next_coupon_date`.

### Accrued interest

```fortran
result = accrued_interest(start_date, end_date, coupon_rate_percent, dcc, &
  redemption_value, coupon_frequency, maturity, next_coupon_year, &
  end_of_month, anniversary_dates)
```

Alias: `accr_int`.

### Clean and dirty price

```fortran
result = dirty_price(clean_price, settlement_date, terms [, schedule])
```

Alias: `dp_value`.

### Price from yield

```fortran
result = bond_price(yield_percent, settlement_date, terms, &
  simple_last_period, calculation_method, supplied_schedule)
```

Alias: `bond_val_price`.

`calculation_method=1` uses the bond's day-count convention for discount-period
coordinates. `calculation_method=0` reproduces the package's alternative ICMA
coordinate calculation. `simple_last_period` defaults to true.

### Yield from price

```fortran
result = bond_yield(clean_price, settlement_date, terms, &
  simple_last_period, precision, calculation_method, supplied_schedule)
```

Alias: `bond_val_yield`. The implementation uses a safeguarded bracketed solve.

## Compatibility module

```fortran
use bondvaluation_compat
```

It exposes Fortran-compatible spellings of the R exports:

- `AccrInt`
- `AnnivDates`
- `DP`
- `BondVal_Price`
- `BondVal_Yield`

It also exposes translated native kernels:

- `leap`, `LDM`, `DaysInMonth`, `DaysInYear`, `DayDiff`, `Date_LDM`
- `sumC`, `FirstMatch`, `LeapDayInside`, `DIST`, `PayCalc`
- `NumToDate`, `CppPrevDate`, `CppSuccDate`
- `NewtonRaphson`, `dm_MyPriceEqn`, `ModDUR`, `CONV`

A dot is not permitted in a Fortran identifier, so R names `BondVal.Price` and
`BondVal.Yield` use underscores.

## Status codes

Schedule codes:

- `0`: success
- `1`: invalid date
- `2`: invalid coupon frequency
- `3`: invalid date ordering
- `4`: schedule unavailable

Pricing additionally uses `10` for invalid settlement range, `11` for no
remaining payment, `12` for an invalid discount base, and `20` if a yield root
cannot be bracketed.
