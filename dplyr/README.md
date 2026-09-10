# dplyr

This package is a restricted modern free-form Fortran translation of **dplyr**
1.2.1, *A Grammar of Data Manipulation*. It operates on the heterogeneous
rectangular types supplied by the sibling `vctrs` and `tibble` translations.

## Implemented foundation

- filtering by a typed logical vector, intrinsic logical array, or expression string;
- selecting, renaming, relocating, mutating, slicing, ordering, and deduplicating;
- row and column binding;
- inner, left, right, full, semi, anti, and cross joins;
- stable grouping, group metadata, counting, and built-in summaries;
- union, union-all, intersection, difference, symmetric difference, and set equality;
- common vector helpers including `between`, `coalesce`, `lag`, `lead`, `near`, and `na_if`.

## Concise filter expressions

The user-facing expression interface supports syntax such as:

```fortran
result = filter(data, "x > 3 and y <= 10")
result = filter(data, "sector == 'energy'")
result = filter(data, "abs(return) < 0.05 and not missing(volume)")
```

Its deliberately small grammar includes parentheses, numeric and character
literals, column names, `+ - * /`, comparisons, `and`/`&`, `or`/`|`,
`not`/`!`, `abs()`, and `missing()`. Both `/=` and `!=` mean not equal. Names
and types are checked at runtime with focused errors.

The statically typed interface uses the same vector operations:

```fortran
result = filter(data, (col(data, "x") > 3.0_dp) .and. &
                      (col(data, "y") <= 10.0_dp))
result = mutate(data, column("z", col(data, "x") + col(data, "y")))
```

Both interfaces use vctrs size-one recycling, numeric promotion, and explicit
missing masks. A missing filtering result discards its row, matching dplyr.

## Grouped summaries

```fortran
grouped = group_by(data, [character(len=6) :: "sector"])
result = summarise(grouped, [ &
   summary("n", "", "n"), &
   summary("average", "return", "mean", na_rm=.true.)])
```

The built-in reductions are `n`, `sum`, `mean`, `min`, and `max` for integer
and double-precision real columns.

## Build and test

```text
fpm build
fpm test
fpm run --example filter_and_summarise
```

The package depends on local `vctrs`, `tibble`, and `tidyselect` packages and
has no external system-library dependency. `select`, `rename`, and `relocate`
accept either exact names or one-based positions returned by tidyselect.

## Scope

This is a portable table-algorithm library, not an implementation of R's data
mask, quosures, S3 dispatch, SQL translation, or lazy database tables. The
translated tidyselect package supplies explicit selectors rather than R's
nonstandard-evaluation expression language. See `API_COVERAGE.md`.

## License

MIT, matching upstream dplyr. See `LICENSE` and `NOTICE.md` for attribution.
