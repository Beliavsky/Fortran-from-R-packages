# Translation coverage

`translation_status: partial`

`frac_functions_translated: 45/284 (0.158)`

The denominator is the 284 functions exported by upstream dplyr 1.2.1. The
numerator counts exported R names with direct Fortran counterparts. It excludes
additional Fortran support procedures and APIs merely re-exported from sibling
packages.

## Direct mappings

| Area | R function | Fortran counterpart |
|---|---|---|
| Rows and columns | `filter` | `filter` |
| Rows and columns | `select` | `select` |
| Rows and columns | `rename` | `rename` |
| Rows and columns | `relocate` | `relocate` |
| Rows and columns | `mutate` | `mutate` |
| Rows and columns | `slice` | `slice` |
| Rows and columns | `slice_head`, `slice_tail` | `slice_head`, `slice_tail` |
| Rows and columns | `arrange` | `arrange` |
| Rows and columns | `distinct` | `distinct` |
| Rows and columns | `pull` | `pull` |
| Binding | `bind_rows`, `bind_cols` | `bind_rows`, `bind_cols` |
| Mutating joins | `inner_join`, `left_join`, `right_join`, `full_join` | Same names |
| Other joins | `semi_join`, `anti_join`, `cross_join` | Same names |
| Grouping | `group_by`, `ungroup`, `n_groups` | Same names |
| Grouping | `group_indices`, `group_keys`, `group_size` | Same names |
| Summaries | `summarise` | `summarise` with `summary_spec` values |
| Summaries | `count` | `count_rows` |
| Sets | `union`, `union_all`, `intersect`, `setdiff` | Same names |
| Sets | `setequal`, `symdiff` | Same names |
| Vector helpers | `between`, `coalesce`, `consecutive_id` | Same names |
| Vector helpers | `first`, `last`, `nth`, `lag`, `lead` | Same names |
| Vector helpers | `n_distinct`, `na_if`, `near` | Same names |

## Additional Fortran API

- `col(data, name)` returns an exact typed column.
- `column(name, expression)` names a computed vector for `mutate()`.
- Arithmetic, comparison, and logical operators are overloaded for `vctr_type`.
- `evaluate_expression()` exposes the filter expression evaluator.
- `grouped_df` and `summary_spec` make grouping and reductions explicit.
- `summary()` constructs one built-in aggregation request.

## Deliberately restricted semantics

- Column selection uses exact names rather than tidyselect expressions.
- Expression strings implement a documented small grammar, not R syntax in general.
- Joins support same-name equality keys. Inequality, rolling, overlap, and differently
  named keys are not yet implemented.
- Paired missing join keys match. Many-to-many matches are allowed and deterministic.
- A conflicting right non-key join column receives `.y`; more complicated name repair
  and caller-selected suffixes are not implemented.
- `arrange()` performs a stable comparison sort and always places missing values last.
- Group summaries provide `n`, `sum`, `mean`, `min`, and `max`, not arbitrary callbacks.
- Table set operations require identical ordered names and types.

## Not translated

- data masking, quosures, tidy evaluation, and scoped/across variants;
- tidyselect helpers and formula interfaces;
- database sources, SQL generation, lazy tables, and distributed backends;
- arbitrary callbacks in `mutate()`, `filter()`, and `summarise()`;
- list columns, rowwise data frames, nesting, and group-map callbacks;
- non-equality join specifications and row-conflict operations;
- deprecated underscore interfaces and compatibility aliases;
- R condition objects, warning registries, and S3 extension hooks.

