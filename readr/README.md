# readr

This package is a restricted modern free-form Fortran translation of
**readr 2.2.0**, *Read Rectangular Text Data*. It provides eager, deterministic
reading and writing of heterogeneous rectangular text files using the sibling
`tibble`, `vctrs`, and `forcats` translations.

## Delimited input

```fortran
type(read_result_type) :: result

result = read_csv("observations.csv")
call stop_for_problems(result)
call print_tibble(result%data)
```

`read_csv`, `read_csv2`, `read_tsv`, `read_delim`, and `read_table` support:

- headers or generated column names;
- explicit named column specifications and skipped columns;
- deterministic logical, integer, real, and character type inference;
- configurable missing tokens;
- quoted delimiters, doubled quotes, and quoted newlines;
- LF and CRLF records, full-line comments, and skipped rows; and
- structured row, column, expected-value, actual-value, and message diagnostics.

Parsing failures become missing values and are retained in the result rather
than silently discarded. `result%ok()`, `problems()`, and
`stop_for_problems()` provide three levels of checking.

## Column parsers

The vector parsers are `parse_logical`, `parse_integer`, `parse_double`,
`parse_number`, `parse_character`, `parse_factor`, `parse_guess`, and
`parse_vector`. Numeric parsing supports configurable decimal and grouping
marks. Factor parsing returns the `factor_type` supplied by the forcats
translation.

## Output

`write_csv`, `write_csv2`, `write_tsv`, and `write_delim` quote fields only
when necessary and double embedded quotes. Their `format_*` counterparts
return `formatted_lines_type`, whose `line` component holds the records.

## Build and test

```text
fpm build
fpm test
fpm run --example read_csv_example
```

Run independent behavior checks against the installed R package with:

```text
Rscript reference/reference_invariants.R
```

## Scope

This is a portable eager parser, not a translation of readr's R interfaces
or its vroom-based lazy, memory-mapped, multithreaded engine. Date/time parsing,
fixed-width files, encodings, compression, URLs, connections, callbacks, RDS,
and clipboard access are deferred. See `API_COVERAGE.md` for exact mappings.

## License

MIT, matching upstream readr. See `LICENSE` and `NOTICE.md` for attribution.
