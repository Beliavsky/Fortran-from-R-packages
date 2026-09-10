# tibble

This package is a modern free-form Fortran translation of the portable data-container
semantics in R package **tibble** 3.3.1, *Simple Data Frames*. It implements a
rectangular heterogeneous table whose integer, double-precision real, logical, and
character columns retain their types. Missingness is represented by a mask on every
column rather than by type-specific sentinel values.

This is deliberately a partial translation. It provides useful static-language
counterparts to construction, validation, exact name lookup, typed extraction,
subsetting, column and row insertion, vector framing, and compact inspection. It does
not emulate R's S3 object system, tidy evaluation, dynamic dots, arbitrary R objects,
or list/matrix columns. See `API_COVERAGE.md` for the exact mapping.

## Build and test

```text
fpm build
fpm test
fpm run --example mixed_tibble
```

The package depends on the sibling `vctrs` translation for typed vector storage,
missingness, slicing, insertion, and type promotion. It has no external dependencies.

## Example

```fortran
use tibble, only: dp, make_column, new_tibble, print_tibble, &
   tibble_column, tibble_type

type(tibble_column) :: columns(3)
type(tibble_type) :: table

columns(1) = make_column("id", [1, 2, 3])
columns(2) = make_column("label", [character(len=5) :: "alpha", "beta", "gamma"])
columns(3) = make_column("score", [2.5_dp, 4.0_dp, 8.5_dp])
table = new_tibble(columns)
call print_tibble(table)
```

`make_column()` also accepts a scalar plus a requested length, reproducing tibble's
size-one recycling rule:

```fortran
table = add_column(table, make_column("selected", .true., table%nrow()))
```

Named access is exact; partial matching is never performed. All columns must have the
same length, and `validate_tibble()` can check a table without terminating the program.

## Design relationship to `r_mod.f90`

The homogeneous `r_tibble_real_t` and `r_tibble_integer_t` implementation in
[`r_mod.f90`](https://github.com/Beliavsky/R-to-Fortran/blob/main/src/r_mod.f90)
informed this package's API and tests. Those
optimized homogeneous representations remain appropriate for generated numerical
code. General heterogeneous vector storage now belongs to the sibling `vctrs` package,
so tibble, future dplyr/tidyr translations, and the transpiler can share one set of
casting and recycling rules without copying the large runtime.

## License

MIT, matching tibble and the borrowed R-to-Fortran design. See `LICENSE` and
`NOTICE.md` for attribution and provenance.
