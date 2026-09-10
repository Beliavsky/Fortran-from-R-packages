# Translation coverage

`translation_status: partial`

`frac_functions_translated: 99/159 (0.623)`

The denominator is the 159 function names exported by upstream lubridate
1.9.5. The numerator counts exported R functions with direct Fortran
counterparts. Types, overloaded operators, and support procedures are not
counted.

## Direct mappings

| Area | R function | Fortran counterpart |
|---|---|---|
| Construction | `make_date`, `make_datetime`, `origin` | Same names |
| Conversion | `as_date`, `as_datetime` | Same names for translated types |
| Date parsing | `ymd`, `ydm`, `mdy`, `myd`, `dmy`, `dym`, `ym`, `my`, `yq` | Same names; numeric fields |
| Date-time parsing | `ymd_h`, `ymd_hm`, `ymd_hms`, `ydm_h`, `ydm_hm`, `ydm_hms`, `mdy_h`, `mdy_hm`, `mdy_hms`, `dmy_h`, `dmy_hm`, `dmy_hms` | Same names; numeric fields |
| Ordered parsing | `parse_date_time` | Same name with one explicit date order |
| Clock parsing | `hm`, `hms`, `ms` | Same names, returning `period_type` |
| Components | `date`, `year`, `month`, `day`, `mday`, `yday`, `qday`, `wday`, `week`, `hour`, `minute`, `second` | Same names |
| Derived calendar fields | `isoweek`, `isoyear`, `epiweek`, `epiyear`, `quarter`, `semester`, `days_in_month`, `leap_year` | Same names |
| Decimal years | `decimal_date`, `date_decimal` | Same names |
| Durations | `duration`, `dseconds`, `dmilliseconds`, `dmicroseconds`, `dnanoseconds`, `dpicoseconds`, `dminutes`, `dhours`, `ddays`, `dweeks`, `dmonths`, `dyears` | Same names |
| Periods | `period`, `seconds`, `milliseconds`, `microseconds`, `nanoseconds`, `picoseconds`, `minutes`, `hours`, `days`, `weeks`, `months`, `years` | Same names |
| Span conversion | `period_to_seconds`, `seconds_to_period`, `time_length` | Same names |
| Intervals | `interval`, `int_start`, `int_end`, `int_length`, `int_flip`, `int_shift`, `int_overlaps`, `int_aligns`, `int_standardize` | Same names |
| Month arithmetic | `add_with_rollback`, `rollback`, `rollforward` | Same names; operators replace `%m+%` and `%m-%` |
| Rounding | `floor_date`, `ceiling_date`, `round_date` | Same names for common clock and calendar units |
| Fixed offsets | `with_tz`, `force_tz` | Same names with integer offsets rather than zone names |
| Clock queries | `now`, `today`, `am`, `pm` | Same names |
| Features | `cyclic_encoding` | Same name |

## Additional Fortran API

- `date_type` and `datetime_type` provide validation and ISO-8601 formatting.
- `duration_type`, `period_type`, and `interval_type` make span semantics
  explicit at compile time.
- `date_to_day_number`, `date_from_day_number`, `datetime_to_epoch`, and
  `datetime_from_epoch` provide deterministic serial conversions.
- Arithmetic and comparison operators cover translated date and span types.
- `within` is the identifier-safe counterpart of R's `%within%` operator.

## Deliberately restricted semantics

- Dates use the proleptic Gregorian calendar, including before 1970.
- Date-times use fixed offsets in minutes east of UTC. Named zones and DST are
  intentionally excluded because they require a changing zone database.
- Parsing extracts numeric components in an explicit order; it does not guess
  locales, month names, truncated forms, or heterogeneous order vectors.
- Calendar month addition clips invalid days to the target month's final day.
- Average duration months and years follow lubridate's 365.25-day year.
- Scalar procedures are elemental where possible, allowing direct array use.

## Not translated

- R class constructors, coercions, predicates, replacement functions, and
  S3/S4 methods tied to `Date`, `POSIXct`, `POSIXlt`, and `difftime`;
- Olson/IANA zones, daylight-saving detection, `force_tzs`, and timeline
  fitting supplied by upstream `timechange`;
- locale-aware names, natural-language parsing, format guessing, stamps,
  `fast_strptime`, and `parse_date_time2`;
- pretty date breaks and R set operations on temporal vectors; and
- formula-style vector recycling, `NA` propagation, and R warning behavior.
