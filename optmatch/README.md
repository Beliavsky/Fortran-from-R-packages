# optmatch-fortran

Modern Fortran translation of the computational core of the R package
`optmatch` 0.10.8, packaged with FPM.

The upstream package implements distance-based optimal bipartite matching for
observational studies, especially optimal full matching through a minimum-cost
flow formulation described by Hansen and Klopfer (2006).

## Scope

This translation concentrates on the numerical and combinatorial work rather
than emulating R's object system.  The public Fortran API accepts ordinary
arrays and derived types.

Implemented computational families include:

- dense and infinity-sparse distance specifications;
- conversion, subsetting, block-diagonal binding, and eligible-edge counting;
- absolute score distances;
- Euclidean distances;
- pooled-covariance Mahalanobis distances;
- rank-based Mahalanobis distances, including the upstream tie adjustment and
  generalized inverse behavior;
- caliper restrictions;
- exact and anti-exact matching restrictions;
- addition/intersection of distance specifications;
- connected matching-subproblem discovery;
- the Hansen-Klopfer full-matching network construction;
- integral minimum-cost flow used by `fullmatch`/`pairmatch`;
- 1:k pair matching and treatment/control omission logic;
- `minControlsCap` and `maxControlsCap` style feasibility searches;
- pooled robust/classical standardization scales;
- exact caliper-size computation and maximum-caliper selection;
- matched-set distances, matched-unit flags, stratum summaries, and effective
  sample size;
- numeric-column NA mean imputation corresponding to the computational numeric
  branch of `fill.NAs`.

The umbrella module is `optmatch`.

## Deliberately not translated

R-specific infrastructure is not useful as Fortran code and is therefore not
ported: S3/S4 classes and dispatch, formulas/model frames, `glm`/`bigglm`/
`survey` integration, data-frame/tibble manipulation, printing, hashing of R
objects, vignettes, and package lifecycle hooks.  Plotting and boxplot code is
omitted as requested.

For model-based matching, pass already-computed fitted scores or a numeric data
matrix to the Fortran distance routines.  This preserves the numerical work
without requiring an R formula interpreter inside the Fortran library.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic_matching
```

No external numerical library is required.  The min-cost-flow solver and the
symmetric pseudoinverse used by the distance calculations are native Fortran.

## Main API mapping

| R `optmatch` computation | Fortran API |
|---|---|
| `fullmatch()` | `fullmatch()` |
| `pairmatch()` / `pair()` | `pairmatch()` |
| internal `fmatch()` | `fmatch_core()` |
| `match_on.numeric()` | `score_distance()` |
| Euclidean `match_on()` distance | `euclidean_distance()` |
| Mahalanobis `match_on()` | `mahalanobis_distance()` |
| rank Mahalanobis | `rank_mahalanobis_distance()` |
| `caliper()` | `caliper_distance()` |
| `exactMatch()` | `exact_match_distance()` |
| `antiExactMatch()` | `anti_exact_match_distance()` |
| distance `+` restriction | `add_distances()` |
| `as.InfinitySparseMatrix()` | `to_sparse()` / `from_sparse()` |
| sparse subsetting C++ helper | `subset_sparse()` |
| `dbind()` | `dbind_distances()` |
| `findSubproblems()` | `find_subproblems()` |
| `num_eligible_matches()` | `num_eligible_matches()` |
| `standardization_scale()` | `standardization_scale()` |
| `effectiveSampleSize()` | `effective_sample_size()` |
| `stratumStructure()` | `summarize_strata()` |
| `matched.distances()` | `matched_distances()` |
| `matched()` / `unmatched()` | `matched_units()` |
| `caliperSize()` | `caliper_size()` |
| `caliperUpperBound()` | `caliper_upper_bound()` |
| `maxCaliper()` | `max_caliper()` |
| `maxControlsCap()` | `max_controls_cap()` |
| `minControlsCap()` | `min_controls_cap()` |
| numeric `fill.column()` | `fill_na_columns()` |
| C++ `digits()` | `integer_digits()` |

`caliper_upper_bound()` intentionally returns the exact eligible-edge count.
That is a sharper valid upper bound than the histogram approximation used by
R and avoids importing plotting/graphics machinery merely to get histogram
counts.

## Distance representation

`type(distance_spec)` stores

- `value(nt, nc)`: nonnegative finite discrepancy values, and
- `allowed(nt, nc)`: whether each treatment-control pairing is eligible.

This is the Fortran analogue of the finite/infinite semantics of
`InfinitySparseMatrix`.  A forbidden edge does not need an IEEE infinity in
`value`; `allowed=.false.` is authoritative.

`type(match_result)` returns set labels separately for treatment and control
units.  Label 0 means unmatched.  `selected(nt,nc)` records the actual
minimum-cost-flow match arcs, while units sharing a positive group label form a
matched set.

## Licensing

The upstream package is MIT licensed.  Its copyright and license are preserved
in `LICENSE`, while the exact CRAN metadata license file is retained in
`original/LICENSE`.  `original/DESCRIPTION` and `original/NAMESPACE` are also
included for provenance.

See `VALIDATION.md` and `TRANSLATION_NOTES.md` for implementation and test
information.
