# Translation coverage

`translation_status: partial`

`frac_functions_translated: 17/70 (0.243)`

The denominator is the 70 functions exported by upstream tidyr 1.3.2. The
numerator counts exported R names with direct Fortran counterparts. It does not
count support procedures or functions re-exported from `tibble` and `vctrs`.

## Direct mappings

| R function | Fortran counterpart | Status |
|---|---|---|
| `pivot_longer()` | `pivot_longer()` | Explicit columns; `fastest` and `slowest` ordering |
| `pivot_wider()` | `pivot_wider()` | One name and one value column; duplicate cells rejected |
| `gather()` | `gather()` | Legacy wrapper around `pivot_longer()` |
| `spread()` | `spread()` | Legacy wrapper around `pivot_wider()` |
| `drop_na()` | `drop_na()` | All columns or an explicit name array |
| `replace_na()` | `replace_na()` | Named scalar replacements with checked casts |
| `fill()` | `fill()` | `down`, `up`, `downup`, and `updown` |
| `separate_wider_delim()` | `separate_wider_delim()` | Character input and exact literal split count |
| `separate()` | `separate()` | Legacy literal-delimiter wrapper |
| `separate_longer_delim()` | `separate_longer_delim()` | Splits one character column into repeated rows |
| `separate_rows()` | `separate_rows()` | Legacy wrapper around `separate_longer_delim()` |
| `unite()` | `unite()` | Formats supported scalar columns and joins them |
| `extract_numeric()` | `extract_numeric()` | Extracts a signed decimal token from text |
| `expand_grid()` | `expand_grid()` | Named `vctr_type` inputs |
| `crossing()` | `crossing()` | Deduplicates and sorts before expansion |
| `full_seq()` | `full_seq()` | Integer and double-precision sequences |
| `uncount()` | `uncount()` | Nonnegative integer weight column |

## Deliberately restricted semantics

- Column selections are exact name arrays rather than tidyselect expressions.
- `pivot_longer()` requires selected columns to have a common `vctrs` type.
- `pivot_wider()` requires a character `names_from` column and rejects repeated
  identifier/name cells instead of accepting an aggregation callback.
- Absent wide cells are represented through the value column's missing mask.
- `separate_wider_delim()` requires exactly the requested number of pieces and
  uses a literal delimiter rather than a regular expression.
- `unite()` formats integer, real, logical, and character columns using compact
  Fortran formatting; exact R display formatting is not promised.
- `expand_grid()` is bounded by the range of the default Fortran integer kind.

## Not translated

- tidyselect helpers and nonstandard evaluation;
- list-column nesting, chopping, packing, hoisting, and rectangling;
- regular-expression extraction and separation;
- pivot specifications, multiple `names_from`/`values_from` columns, and
  callback-based aggregation;
- deprecated underscore interfaces;
- arbitrary R classes and S3 methods.
