# Translation coverage

`translation_status: partial`

`frac_functions_translated: 20/184 (0.109)`

The denominator is the 184 names exported by purrr 1.2.2. The numerator counts
exported R functions with direct typed Fortran counterparts.

| R function family | Fortran counterpart |
|---|---|
| `map`, `map_dbl`, `map_int`, `modify` | Generic and explicitly typed real/integer maps |
| `map2`, `map2_dbl`, `map2_int` | Conformable binary maps |
| `reduce` | Typed left folds with optional initial value |
| `accumulate` | Typed intermediate left-fold results |
| `keep`, `discard` | Predicate-based typed filtering |
| `some`, `every`, `none` | Predicate quantifiers |
| `detect`, `detect_index` | First matching value or position |
| `head_while`, `tail_while` | Contiguous predicate-based slices |
| `walk` | Side-effecting typed traversal |
| `negate` | Direct predicate negation for one real value |

The callback contracts are public abstract interfaces. Pure operations reject
impure callbacks at compile time.

## Not translated

- heterogeneous R lists, quosures, formulas, dynamic dots, and environments;
- character/raw/data-frame variants and arbitrary output-type simplification;
- `pmap`, depth traversal, list mutation, plucking, splicing, and transposition;
- error capture and retry tools such as `safely`, `possibly`, and `insistently`;
- parallel and asynchronous mapping; and
- R type predicates that have no role in statically typed Fortran.
