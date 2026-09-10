# forcats

This package is a restricted modern free-form Fortran translation of
**forcats 1.0.1**, *Tools for Working with Categorical Variables (Factors)*.
It provides a type-safe representation of categorical observations and the
portable deterministic algorithms at the center of the R package.

## Factor representation

`factor_type` stores:

- one-based integer observation codes;
- an ordered array of character level labels;
- an explicit missing-value mask; and
- whether the factor is ordinal.

Constructors validate code ranges, level uniqueness, and mask conformance.
Character conversion uses first-appearance order, so it is deterministic
across locales. Integer and real conversion orders levels numerically.

```fortran
type(factor_type) :: sector, compact

sector = fct([character(len=10) :: "technology", "finance", &
   "technology", "energy"])
compact = fct_lump_min(sector, 2.0_dp, other_level="other")
```

## Implemented operations

- construct and decode factors;
- reorder levels manually, by appearance, frequency, numeric value, cyclic
  shift, or aligned numeric summaries;
- recode, collapse, retain, expand, and drop levels;
- convert between missing observations and an explicit missing level;
- concatenate, unify, cross, match, and enumerate factors;
- count levels and compute proportions; and
- lump levels by count, proportion, rank, or collective low frequency.

All library procedures are pure. Random shuffling and anonymization are not
included because they require mutable random-number state; callers can create
a permutation with their chosen generator and pass it to `lvls_reorder()`.

## Build and test

```text
fpm build
fpm test
fpm run --example factor_workflow
```

Run independent behavior checks against the installed R package with:

```text
Rscript reference/reference_invariants.R
```

## Scope

Fortran uses explicit arrays in place of R dynamic dots and accepts named
built-in summaries in place of arbitrary R callbacks. `fct_cross()` currently
accepts two factors. The package does not reproduce R's S3 dispatch,
attributes, tidy evaluation, condition objects, or pipe operator. See
`API_COVERAGE.md` for exact mappings.

## License

MIT, matching upstream forcats. See `LICENSE` and `NOTICE.md` for attribution.
