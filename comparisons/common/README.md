# Shared comparison data support

This FPM package provides reusable dated asset-price data support for the
deterministic R-versus-Fortran comparison suites.

`asset_price_data` stores only observed dates, symbols, and prices. Explicit
`simple_returns` and `log_returns` transformations produce a separate
`asset_return_data` value whose dates are the ending dates of the return
intervals. Data loading and return construction should occur outside timed
sections of comparisons.

The human-readable CSV file is the canonical source. `convert-prices` creates a
versioned binary fixture with an explicitly little-endian layout: a signature,
32-bit dimensions, fixed-width symbols, 32-bit year/month/day components, and
IEEE 64-bit prices. Both the Fortran reader and `read_asset_prices_binary.R`
read this format without relying on compiler-specific derived-type layout.

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

Benchmark CSV reading and return construction with:

```text
fpm run --profile release --example benchmark-read
```

Regenerate the binary fixture from the CSV with:

```text
fpm run --profile release --example convert-prices
```

The tests read `../../asset_class_etf_prices.csv`, verify its dimensions, date
range, symbols, and first simple and logarithmic returns, and require the binary
reader to reproduce every CSV date, symbol, and price exactly.
