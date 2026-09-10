# tidyr

This package is a restricted modern free-form Fortran translation of **tidyr**
1.3.2, *Tidy Messy Data*. It provides deterministic transformations of the
heterogeneous rectangular data types implemented by the sibling `vctrs` and
`tibble` translations.

## Implemented foundation

- `pivot_longer()` and `pivot_wider()` for explicit columns;
- legacy `gather()` and `spread()` wrappers;
- `drop_na()`, `replace_na()`, and four-direction `fill()`;
- exact literal-delimiter wider and longer separation, plus `unite()`;
- `expand_grid()` and sorted, deduplicated `crossing()`;
- `full_seq()` for aligned integer and real sequences; and
- integer-weight `uncount()`, with optional within-row identifiers.

Column selections use exact character names. This is intentional: Fortran has
no direct counterpart to tidyselect's R expression language. Missing values use
the explicit masks in `vctr_type`, so valid integer and character values never
double as sentinels.

## Build and test

```text
fpm build
fpm test
fpm run --example reshape_measurements
```

The package uses the local `vctrs`, `tibble`, `tidyselect`, and `stringr`
translations. `pivot_longer` accepts either exact names or positions returned
by tidyselect. Literal wider splitting reuses stringr's tested splitter.

## Example

```fortran
use tidyr, only: dp, make_column, new_tibble, pivot_longer, &
   tibble_column, tibble_type

type(tibble_column) :: columns(3)
type(tibble_type) :: long, wide

columns(1) = make_column("station", [1, 2])
columns(2) = make_column("morning", [12.5_dp, 14.0_dp])
columns(3) = make_column("evening", [10.0_dp, 11.5_dp])
wide = new_tibble(columns)
long = pivot_longer(wide, [character(len=7) :: "morning", "evening"])
```

## Scope

This is a portable data-transformation library, not an implementation of R's
evaluation model. It excludes R quosures and data-mask expressions, list-columns,
nested data frames, regular expressions, and arbitrary callback functions. See
`API_COVERAGE.md` for exact mappings and restrictions.

## License

MIT, matching upstream tidyr. See `LICENSE` and `NOTICE.md` for attribution.
