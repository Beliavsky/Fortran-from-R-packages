# API

All public names are available from module `creditr`.

## Types

- `date_t`: Gregorian date with `year`, `month`, `day`, and `serial()`.
- `cds_dates_t`: trade, step-in, value, start, first-coupon, penultimate-coupon, end, backstop, and base dates.
- `conventions_t`: day counts, frequencies, calendars, and bad-day convention.
- `rate_quote_t`: expiry in months, instrument type (`M` or `S`), and decimal market rate.
- `zero_curve_t`: base date, node dates, log discount factors, `discount`, and `forward_discount`.
- `cds_contract_t`: CDS terms.
- `cds_result_t`: clean/dirty values, price, probability, hazard, and sensitivities.

## Date and convention routines

```fortran
make_date(year, month, day)
add_dates(trade_date, currency, dates, tenor_years=..., maturity=..., status=...)
add_conventions(currency, conventions, status)
```

Exactly one of `tenor_years` or `maturity` must be supplied to `add_dates`.

## Curve routines

```fortran
read_rate_quotes_csv(filename, quotes, status)
build_zero_curve(base_date, quotes, conventions, curve, status)
df = curve%discount(date)
fwd = curve%forward_discount(start_date, end_date)
```

Rates are decimal annual rates. Expiries in CSV use labels such as `1M`, `6M`, `2Y`, and `30Y`.

## CDS valuation

```fortran
result = price_cds(contract, curve, status, quotes)
result = cds(contract, curve, status, quotes)
```

The optional `quotes` argument is used for exact market-quote interest-rate DV01.

```fortran
upfront = spread_to_upfront(curve, dates, spread_bps, coupon_bps, &
  recovery, notional, is_price_clean, status)
spread = upfront_to_spread(curve, dates, upfront, coupon_bps, &
  recovery, notional, is_price_clean, status)
```

## Risk measures

```fortran
value = spread_dv01(curve, dates, spread_bps, coupon_bps, recovery, notional)
value = cs10(curve, dates, spread_bps, coupon_bps, recovery, notional)
value = rec_risk_01(curve, dates, spread_bps, coupon_bps, recovery, notional)
value = ir_dv01(curve, dates, spread_bps, coupon_bps, recovery, notional)
```

The standalone `ir_dv01` applies a parallel zero-rate bump. `price_cds(..., quotes=...)` rebuilds from bumped input quotes as the R package does.

## Simple formulas

```fortran
pd = spread_to_pd(spread_bps, recovery, time_years)
spread_bps = pd_to_spread(pd, recovery, time_years)
recovery_percent = implied_rr(pd, spread_bps, tenor_years)
value = pv01(principal, notional, spread_bps, coupon_bps)
```
