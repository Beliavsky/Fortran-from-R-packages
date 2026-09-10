# Translation coverage

`translation_status: partial`

`frac_functions_translated: 12/45 (0.267)`

The denominator is the 45 functions exported in the upstream tibble 3.3.1
`NAMESPACE`. Multiple upstream entry points can map to one typed Fortran API; the
fraction measures exported-name coverage, not line coverage or semantic completeness.

## Mapped upstream functions

| R function or method | Fortran counterpart | Scope |
|---|---|---|
| `tibble()`, `new_tibble()` | `new_tibble()` | Construct from explicitly typed columns; optional row count for zero-column tables |
| `as_tibble()` | `as_tibble()` | Integer, real, logical, or character rank-2 arrays |
| `validate_tibble()` | `validate_tibble()` | Rectangular shape, storage tag, missing mask, and name invariants |
| `repair_names()` | `repair_names()`; `new_tibble(..., name_repair=)` | `check_unique`, `unique`, and `minimal` modes |
| `has_name()` | `table%has_name()` | Exact matching only |
| `[` | `slice_rows()`, `filter_rows()`, `slice_columns()`, `select_columns()` | Explicit typed indices or masks; preserves table shape |
| `[[` | `get_integer()`, `get_real()`, `get_logical()`, `get_character()` | Exact typed named extraction |
| `[[<-` | `replace_column()` | Existing named column; its type may change |
| `add_column()` | `add_column()` | Explicit insertion position; duplicate names rejected |
| `add_row()`, `add_case()` | `add_rows()` | Inserts a name-compatible table, with vctrs type promotion |
| `enframe()` | generic `enframe()` | Integer, real, logical, and character vectors |
| `deframe()` | `deframe_integer()`, `deframe_real()`, `deframe_logical()`, `deframe_character()` | One- or two-column tables with typed results |
| `glimpse()` | `glimpse()` | Compact transposed preview |
| `print.tbl()` / `trunc_mat()` concept | `print_tibble()` | Compact values and tibble-style type abbreviations |

The 12/45 numerator counts exported upstream names, not the additional S3 operator
methods shown in the table.

## Additional Fortran API

- `make_column()` constructs vector columns or recycles scalar values.
- `get_missing()` extracts a named column's missing-value mask independently of type.
- `drop_columns()` removes exact names.
- `tibble_column%type_name()` returns `int`, `dbl`, `lgl`, or `chr`.
- Every column carries a logical missing-value mask, including character and integer
  columns where IEEE NaN is unavailable or inappropriate.

## Deliberately deferred

- R S3 class inheritance and method dispatch
- tidy evaluation, quosures, dynamic dots, and expression-dependent column creation
- arbitrary list, nested tibble, matrix, array, and user-defined R-object columns
- `tribble()`'s variadic formula syntax and `lst()` expression naming
- R row-name conversion helpers; this Fortran type intentionally has no row names
- `view()` and other IDE/interactive behavior
- pillar's full width-, color-, locale-, and terminal-aware formatting
- compatibility aliases for superseded APIs such as `data_frame()` and `as.tibble()`
- re-exported functions owned by magrittr, pillar, rlang, or vctrs

These omissions are language-runtime features rather than missing numerical kernels.
