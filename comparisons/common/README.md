# Shared comparison data support

This FPM package provides reusable dated asset-price data support for the
deterministic R-versus-Fortran comparison suites.

`asset_price_data` stores only observed dates, symbols, and prices. Explicit
`simple_returns` and `log_returns` transformations produce a separate
`asset_return_data` value whose dates are the ending dates of the return
intervals. Data loading and return construction should occur outside timed
sections of comparisons.

The date implementation supports validated ISO dates, date arithmetic and
comparison, ISO weekdays, day-of-year values, quarters, and month-end tests.

Other comparison packages can use this package with:

```toml
[dependencies]
comparison_data = { path = "../common" }
```

Run its tests from this directory:

```text
fpm test
```

The tests read `../../asset_class_etf_prices.csv` and verify its dimensions,
date range, symbols, and first simple and logarithmic returns.
