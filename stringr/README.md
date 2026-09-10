# stringr

This package is a restricted modern free-form Fortran translation of
**stringr 1.6.0**, *Simple, Consistent Wrappers for Common String Operations*.
It provides deterministic fixed-string matching and common transformations
without requiring `stringi`, ICU, or a regular-expression system library.

Most scalar operations are pure elemental functions, so the same call works
on a scalar or a conformable array:

```fortran
character(len=12) :: names(3)

names = str_to_snake(names)
print *, str_detect(names, 'rate')
```

Implemented areas include detection, counting, location, extraction,
substrings, ASCII case conversion, trimming, whitespace squishing, padding,
truncation, replacement, SQL-style wildcard matching, word extraction,
joining, ordering, sorting, ranking, uniqueness, subsets, and literal splits.

## Build and test

```text
fpm build
fpm test
fpm run --example basic
```

## Scope

Patterns are literal fixed strings except for `str_like()`, which implements
SQL `%` and `_` wildcards. Case conversion and character classification use
ASCII rules. Bytes outside ASCII are preserved but are not Unicode-normalized,
case-folded, collated, or treated as grapheme clusters. See `API_COVERAGE.md`.

## License

MIT, matching upstream stringr. See `LICENSE`, `NOTICE.md`, and `upstream/`.
