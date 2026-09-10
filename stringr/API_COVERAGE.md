# Translation coverage

`translation_status: partial`

`frac_functions_translated: 42/63 (0.667)`

The denominator is the 63 names exported by stringr 1.6.0. The numerator
counts exported R functions with direct Fortran counterparts.

## Direct mappings

| Area | R functions | Fortran status |
|---|---|---|
| Size | `str_length`, `str_width` | Trimmed byte length and width |
| Fixed matching | `str_detect`, `str_starts`, `str_ends`, `str_count`, `str_equal` | Pure elemental literal matching |
| Location/extraction | `str_locate`, `str_extract` | First literal match |
| Substrings | `str_sub` | Positive and negative inclusive indices |
| Replacement | `str_replace`, `str_replace_all`, `str_remove`, `str_remove_all` | Literal patterns |
| Case | `str_to_lower`, `str_to_upper`, `str_to_title`, `str_to_sentence` | ASCII case rules |
| Identifier case | `str_to_snake`, `str_to_kebab`, `str_to_camel` | Common ASCII separators |
| Whitespace | `str_trim`, `str_squish` | ASCII whitespace |
| Shape | `str_dup`, `str_pad`, `str_trunc` | Same names |
| Joining | `str_c`, `str_flatten`, `str_flatten_comma` | Explicit string arrays |
| Splitting | `str_split`, `str_split_1`, `str_split_fixed`, `str_split_i` | Literal delimiters |
| Vector operations | `str_order`, `str_rank`, `str_sort`, `str_unique`, `str_subset`, `str_which` | Stable bytewise ordering |
| Wildcards | `str_like` | SQL `%` and `_` patterns |
| Utilities | `str_escape`, `word` | Same names with restricted semantics |

## Additional Fortran API

- `split_result_type` stores a rectangular split result and per-row counts.
- `location_type` stores one-based inclusive start and end positions.
- Scalar fixed-string operations are elemental where their result shape permits.

## Not translated

- ICU collation, locale-aware case conversion, normalization, and conversion;
- regular-expression patterns, capture groups, match matrices, and callbacks;
- boundary iterators, locale collation objects, and pattern wrapper objects;
- glue interpolation, HTML views, arbitrary encodings, and text wrapping; and
- R missing-value recycling and warning behavior.
