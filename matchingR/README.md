# matchingR-fortran 0.1.0

Modern Fortran 2018/FPM translation of the computational algorithms in
`matchingR` 2.0.0.

The library is native Fortran. It does not require R, Rcpp, Armadillo, BLAS, or
another matching library.

## Implemented algorithms

- Gale-Shapley deferred acceptance for one-to-one stable matching.
- College admissions / many-to-one deferred acceptance with heterogeneous
  college capacities, in both student-optimal and college-optimal directions.
- Irving's stable-roommates algorithm, including the dummy participant used by
  matchingR for odd market sizes.
- Top trading cycles for one-sided indivisible-good exchange.
- Stability checks for two-sided matching, stable roommates, and TTC output.
- `sortIndex`, `sortIndexOneSided`, `rankIndex`, and preference-completeness
  functionality in Fortran form.
- Cardinal-utility and ordinal-preference entry points.

## Indexing

The public Fortran API uses normal Fortran **1-based IDs**. The integer value
`0` represents an unmatched participant or vacant college slot. Ordinal input
routines accept either 1-based or matchingR/Rcpp-style 0-based preferences and
normalize them internally.

## Main API

```fortran
use matchingr

type(marriage_result_t) :: m
type(college_result_t)  :: c
type(roommate_result_t) :: r

m = marriage_market(proposer_utils, reviewer_utils)
m = marriage_market_preferences(proposer_pref, reviewer_pref)

c = college_admissions(student_utils, college_utils, slots, &
                       student_optimal=.true.)

r = stable_roommates(utils)
r = stable_roommates_preferences(pref)

matching = top_trading_cycles(utils)
matching = top_trading_cycles_preferences(pref)
```

For `college_result_t`, `matched_colleges(college,slot)` is zero for a vacant
slot and `matched_students(student)` is zero for an unmatched student.

## Build

```text
fpm build
fpm test
fpm run --example marriage_example
```

No external libraries are required.

## Validation

Seven regression executables cover exact examples from matchingR, both
student- and college-optimal deferred acceptance, 0/1-based preference input,
Irving's saved even/odd examples, TTC, utility/rank helpers, and randomized
stability/property checks. See `VALIDATION.md`.

## Compatibility note: roommate stability checker

The original C++ `cpp_wrapper_irving_check_stability()` contains two apparent
source bugs: it can count the current partner as a blocking candidate, and the
second participant's scan compares against that participant's own ID instead of
the other member of the candidate pair. The R wrapper also passes differently
indexed preference/matching values. The Fortran function `roommate_stable()`
uses the mathematically correct blocking-pair test by default. Setting
`legacy_cpp_bug=.true.` reproduces the two low-level C++ scan-order/comparison
quirks after translating IDs to Fortran indexing. The Irving solver itself is a
translation of the matching algorithm, not of this checker bug.

## Scope

R list construction, warnings/errors, Rcpp registration, roxygen/S3 plumbing,
and vignette presentation are omitted. They are interface code rather than
matching algorithms.

## License

GPL-2.0-or-later, matching the upstream `License: GPL (>= 2)` declaration.
The complete attached matchingR source tree is retained under
`original/matchingR-master/` for provenance.
