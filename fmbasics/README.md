# fmbasics-fortran

A modern Fortran translation of the computational core of the R package
`fmbasics` 0.3.99, packaged for the Fortran Package Manager (FPM).

The library provides typed financial-market building blocks:

- civil dates, day-count fractions, periods, business-day calendars and rolls;
- currencies, currency pairs, FX value dates and benchmark rate indices;
- interest rates, discount factors and compounding conversions;
- zero curves and constant, linear, cubic and log-discount interpolation;
- CDS specifications, survival probabilities, hazard rates and credit curves;
- volatility quotes and time-variance volatility-surface interpolation;
- single- and multicurrency money values and dated cash flows.

The implementation is self-contained and has no external FPM dependencies.
R printing, tibble conversion, S3 dispatch and external holiday-database
infrastructure are not translated.

## Build

```text
fpm build
fpm test
fpm run
```

The default executable is `demo_fmbasics`. Examples may be run with:

```text
fpm run --example rates_curve_example
fpm run --example conventions_money_example
fpm run --example credit_vol_example
```

GNU Fortran users can also run:

```text
./scripts/test_gfortran.sh
```

or on Windows:

```text
scripts\test_gfortran.bat
```

## Basic use

```fortran
use fmbasics

type(interest_rate_t) :: rate
type(discount_factor_t) :: df
integer :: d1, d2

d1 = make_date(2020, 1, 1)
d2 = make_date(2025, 1, 1)
rate = interest_rate(0.04_dp, 2.0_dp, 'act/365')
df = as_discount_factor(rate, d1, d2)
```

Dates are stored as integer days relative to 1970-01-01. The helpers
`make_date`, `date_from_yyyymmdd`, `date_to_yyyymmdd` and `date_string`
provide conversion to and from civil dates.

## Compatibility notes

The original package delegates calendar logic to `fmdates` and CDS curve
construction to `credule`. This port supplies deterministic internal
implementations. Weekend rules and a selected set of common holidays are
included, but the calendar database is not a complete replacement for an
actively maintained market-holiday service. The CDS bootstrap assumes
piecewise-constant hazard rates and quarterly premium payments.

See `PORTING.md` and `TRANSLATION_COVERAGE.md` for detailed differences.

## License

The translated project is distributed under GPL-2.0-only, matching the
original package metadata. Original R sources, tests and metadata are retained
under `original/` for provenance and license compliance.
