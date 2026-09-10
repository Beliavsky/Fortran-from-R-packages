# lubridate

This package is a restricted modern free-form Fortran translation of
**lubridate 1.9.5**, *Make Dealing with Dates a Little Easier*. It provides
strongly typed dates, fixed-offset date-times, durations, calendar periods,
intervals, numeric parsing, accessors, calendar arithmetic, and rounding.

```fortran
use lubridate

type(date_type) :: due
type(datetime_type) :: instant

due = ymd("2024-01-31") + months(1) ! 2024-02-29
instant = ymd_hms("2024-06-19 12:34:56")
instant = floor_date(instant, "hour")
```

Date and date-time values are validated proleptic-Gregorian components.
Date-times carry a numeric fixed offset in minutes east of UTC, so conversions
are deterministic and do not depend on a platform time-zone database.
Durations represent exact seconds; periods retain calendar components. This
distinction makes `ddays(1.0_dp)` an exact 86,400 seconds while `days(1)` is a
calendar day.

## Build and test

```text
fpm build
fpm test
fpm run --example basic
```

The tests cover leap and century years, dates before 1970, ISO week boundaries,
weekday conventions, parsing orders, end-of-month rollback, exact durations,
periods, intervals, rounding, and fixed-offset conversion. Independent
reference invariants can be checked with:

```text
Rscript reference/reference_invariants.R
```

## Scope

The translation targets portable deterministic computational behavior. It does
not reproduce R's `POSIXct`/`POSIXlt` classes, vector recycling, S3/S4 dispatch,
locale-dependent names, natural-language parsing, Olson/IANA named zones,
daylight-saving transitions, or `timechange` integration. Parsers currently
accept numeric component strings. See `API_COVERAGE.md` for exact mappings.

## License

MIT, matching upstream lubridate. See `LICENSE`, `NOTICE.md`, and the retained
upstream metadata under `upstream/`.
