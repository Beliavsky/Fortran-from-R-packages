# Translation coverage

`translation_status: partial`

`frac_functions_translated: 46/103 (0.447)`

The denominator is the 103 names exported by upstream readr 2.2.0. The
numerator counts exported R functions with direct Fortran counterparts.
Additional result types and support procedures are not counted.

## Direct mappings

| Area | R function | Fortran counterpart |
|---|---|---|
| Collectors | `col_logical`, `col_integer`, `col_double`, `col_number`, `col_character`, `col_factor`, `col_guess`, `col_skip` | Same names |
| Specifications | `cols` | Same name with parallel name and collector arrays |
| Locales | `locale`, `default_locale` | Same names for decimal/grouping marks and whitespace |
| Tokenization | `tokenize`, `count_fields` | Same names |
| Vector parsing | `parse_logical`, `parse_integer`, `parse_double`, `parse_number` | Same names |
| Vector parsing | `parse_character`, `parse_factor`, `parse_guess`, `parse_vector` | Same names |
| Delimited input | `read_delim`, `read_csv`, `read_csv2`, `read_tsv`, `read_table` | Same names |
| Whole-file input | `read_file`, `read_lines` | Same names |
| Specifications | `spec`, `spec_delim`, `spec_csv`, `spec_csv2`, `spec_tsv`, `spec_table` | Same names |
| Diagnostics | `problems`, `stop_for_problems` | Same names |
| Formatting | `format_delim`, `format_csv`, `format_csv2`, `format_tsv` | Same names |
| Output | `write_delim`, `write_csv`, `write_csv2`, `write_tsv` | Same names |
| Whole-file output | `write_file`, `write_lines` | Same names |

## Additional Fortran API

- `read_result_type` owns a tibble, its effective specification, and problems.
- `parse_result_type` and `factor_parse_result_type` pair values with problems.
- `parse_problem_type` provides stable structured diagnostics.
- `token_table_type` exposes tokens, quoted-field flags, and field counts.
- `formatted_lines_type` owns preformatted output records.

## Deliberately restricted semantics

- Readers are eager and single-threaded.
- A factor collector is represented as a character tibble column because the
  current `vctrs` storage has no categorical column type; `parse_factor()`
  returns a full forcats `factor_type` when level metadata is required.
- Column specifications use explicit parallel arrays rather than R dynamic dots.
- Empty strings and `NA` are missing by default.
- Name repair is deterministic and uses the translated tibble implementation.
- `format_*` returns an array of records rather than one R character scalar.

## Not translated

- date, time, and datetime collectors and locale-specific date names;
- fixed-width and web-log readers;
- vroom lazy parsing, memory mapping, threading, and chunk callbacks;
- R connections, URLs, compression, clipboard, raw vectors, and RDS;
- encoding detection and transcoding;
- readr editions, progress display, global options, and R condition classes;
- tokenizer configuration objects and R S3 formatting methods.
