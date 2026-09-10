# Translation coverage

`translation_status: partial`

`frac_functions_translated: 32/200 (0.160)`

The denominator is the 200 functions exported by upstream vctrs 0.7.3. The numerator
counts exported names with direct Fortran counterparts. It does not count additional
Fortran support routines or upstream S3 methods.

## Directly mapped exports

| Area | Upstream exports with Fortran counterparts |
|---|---|
| Construction and validation | `new_vctr`, `new_data_frame`, `vec_is`, `vec_size` |
| Missingness | `vec_detect_missing`, `vec_any_missing` |
| Types and casting | `vec_ptype2`, `vec_ptype_common`, `vec_cast` |
| Sizes and recycling | `vec_size_common`, `vec_recycle`, `vec_recycle_common` |
| Combining and slicing | `vec_c`, `vec_slice`, `vec_rep`, `vec_rep_each` |
| Comparison and ordering | `vec_equal`, `vec_compare`, `vec_order`, `vec_sort` |
| Matching and sets | `vec_match`, `vec_in`, `vec_unique`, `vec_unique_loc`, `vec_unique_count` |
| Duplicates and grouping | `vec_duplicate_any`, `vec_duplicate_detect`, `vec_group_id`, `vec_identify_runs` |
| Conditional selection | `vec_if_else` |
| Data frames | `vec_rbind`, `vec_cbind` |

## Additional Fortran API

- `vctr_type` is the explicit tagged-vector representation.
- `vctrs_data_frame` is a rectangular collection of named vectors.
- `validate_data_frame()` returns a status and optional diagnostic message.
- `vec_insert()` inserts a compatible vector at a specified position.
- `vec_init()` allocates same-type output storage with an explicit initial mask.
- `vec_copy_element()` copies one typed value and its missingness between vectors.
- `dp` and the four public type-code constants make interfaces explicit.

## Deliberately restricted semantics

- Only logical, default integer, `real(dp)`, and character vectors are represented.
- Common numeric types follow logical to integer to real promotion.
- Character and numeric vectors have no implicit common type.
- Narrowing casts must be lossless; there is no `allow_lossy_cast` escape hatch.
- `vec_rbind()` currently requires identical column names and order. It promotes
  compatible column types but does not construct the union of disjoint schemas.
- `vec_cbind()` requires globally unique names and recycles only one-row frames.
- Sorting is stable and deterministic but uses a simple comparison sort rather than
  upstream's optimized radix implementation.

## Deliberately deferred

- S3 method dispatch, class proxy/restore hooks, and arbitrary user-defined classes
- list, record, factor, ordered-factor, raw, complex, date, datetime, duration, and
  integer64 vector classes
- multidimensional arrays and shaped vector broadcasting
- R expressions, environments, dots, quosures, and condition/error classes
- ALTREP and R/C registration machinery
- locale-aware character collation
- optimized hashing and radix sorting
- data-frame schema union and sophisticated relationship diagnostics for joins

These areas can be added incrementally when a translated downstream package requires
them. They should not be inferred from the package name.
