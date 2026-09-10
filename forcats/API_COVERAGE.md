# Translation coverage

`translation_status: partial`

`frac_functions_translated: 35/39 (0.897)`

The denominator is the 39 names exported by upstream forcats 1.0.1. The
numerator counts exported R functions with a direct Fortran counterpart.
Additional constructors, types, and support procedures are not counted.

## Direct mappings

| Area | R function | Fortran counterpart |
|---|---|---|
| Construction | `fct`, `as_factor` | Generic procedures with character, integer, real, and logical inputs |
| Low-level levels | `lvls_reorder`, `lvls_revalue`, `lvls_expand`, `lvls_union` | Same names |
| Manual ordering | `fct_relevel`, `fct_rev`, `fct_shift` | Same names |
| Data ordering | `fct_inorder`, `fct_infreq`, `fct_inseq` | Same names |
| Summary ordering | `fct_reorder`, `fct_reorder2`, `first2`, `last2` | Same names with built-in summaries |
| Relabeling | `fct_recode`, `fct_collapse`, `fct_other` | Same names with explicit mapping arrays |
| Level sets | `fct_expand`, `fct_drop` | Same names |
| Missing values | `fct_na_value_to_level`, `fct_na_level_to_value` | Same names using explicit masks |
| Deprecated alias | `fct_explicit_na` | Same name, forwarding to `fct_na_value_to_level` |
| Combining | `fct_c`, `fct_unify`, `fct_cross` | Same names; crossing currently accepts two factors |
| Inspection | `fct_count`, `fct_match`, `fct_unique` | Same names |
| Lumping | `fct_lump`, `fct_lump_min`, `fct_lump_prop`, `fct_lump_n`, `fct_lump_lowfreq` | Same names |

## Additional Fortran API

- `factor_type` is the validated categorical container.
- `factor_count_type` holds factor labels, counts, and optional proportions.
- `new_factor()` constructs directly from codes, levels, and a missing mask.
- Type-bound `size()`, `nlevels()`, `values()`, and `valid()` inspect factors.

## Deliberately restricted semantics

- Mapping arrays replace dynamic dots in recoding and collapsing operations.
- `fct_reorder()` supports `mean`, `median`, `min`, and `max` rather than an
  arbitrary callback.
- `fct_reorder2()` provides the upstream first/last terminal-value behavior.
- `fct_lump_n()` preserves all frequency ties at its boundary; configurable
  R rank tie methods are not reproduced.
- Missing values use a separate logical mask rather than an `NA` integer code
  or an `NA_character_` level.
- The compact count-table type replaces an R tibble return value.

## Not translated

- `%>%`, which is R syntax rather than a numerical operation;
- `fct_anon` and `fct_shuffle`, whose implicit global RNG conflicts with the
  pure deterministic API;
- `fct_relabel`, which accepts an arbitrary R function;
