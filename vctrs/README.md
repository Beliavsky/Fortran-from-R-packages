# vctrs

This package is a restricted modern free-form Fortran translation of **vctrs**
0.7.3, *Vector Helpers*. It provides explicit type- and size-stable operations for
logical, integer, double-precision real, and character vectors.

Every `vctr_type` stores exactly one typed value array and a separate logical
missing-value mask. This avoids overloading valid integer or character values as
sentinels and lets downstream packages use the same missingness rules for every type.

## Implemented foundation

- vector construction, validation, size, type names, and missingness;
- common-type determination and checked lossless casting;
- common-size determination and size-one recycling;
- slicing, concatenation, insertion, repetition, and stable sorting;
- equality, three-way comparison, matching, and membership;
- unique values, duplicate detection, grouping, and run identification;
- type-stable conditional selection;
- rectangular data-frame construction, row binding, and column binding.

The promotion hierarchy is:

```text
logical -> integer -> real
```

Character vectors combine only with character vectors. Narrowing casts are accepted
only when every nonmissing value can be represented without loss.

## Build and test

```text
fpm build
fpm test
fpm run --example stable_vectors
```

The package has no dependencies. The sibling `tibble` translation depends on it and
re-exports `vctr_type` and `vctrs_data_frame` under tibble-facing type names.

## Example

```fortran
use vctrs, only: dp, new_vctr, vctr_type, vec_c

type(vctr_type) :: inputs(2), combined

inputs(1) = new_vctr("value", [1, 2])
inputs(2) = new_vctr("value", [3.5_dp])
combined = vec_c(inputs)
```

`combined` is a real vector containing `1.0`, `2.0`, and `3.5` because the common
type of integer and real is real.

## Scope

The package implements portable vector semantics, not R's runtime. S3 dispatch,
arbitrary R classes, language objects, environments, ALTREP, and R condition classes
are deliberately excluded. See `API_COVERAGE.md` for the exact mapping.

## License

MIT, matching upstream vctrs. See `LICENSE` and `NOTICE.md` for attribution.
