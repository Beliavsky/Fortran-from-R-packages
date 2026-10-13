# exactRankTests for Fortran

Modern free-form Fortran translation of exactRankTests **0.8-37**, maintained
upstream by Torsten Hothorn and Kurt Hornik. Provides exact conditional
permutation distributions with ties, rank scoring, Wilcoxon, Ansari-Bradley,
and general permutation tests. See [NOTICE.md](NOTICE.md) for provenance and
[LICENSE](LICENSE) for GPL-2.0-or-later terms.

## Build and run

Place this directory at the top level of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
alongside the existing `rfortran-core` package, then run from this directory:

```sh
fpm build
fpm test
fpm run --example --all
```

Requires FPM and a modern Fortran compiler. The manifest uses the sibling
`rfortran-core` path dependency for `dp`, ordering, average ranks and normal
probabilities/quantiles. No dependency sources are copied, no system BLAS or
LAPACK is linked, and no operating-system-specific APIs are used. Windows with
gfortran/FPM uses the same commands and paths; validation here ran on Linux,
not a Windows host. Never use `-ffast-math`, `-Ofast`, or
`-ffinite-math-only`: error handling depends on IEEE NaNs/infinities.

## Main API

Import `exact_rank_tests` for the public facade; it re-exports `dp`, status
constants, distribution/scoring APIs, and the test-result type.

```fortran
use exact_rank_tests
implicit none
type(rank_test_result) :: result
call wilcox_exact([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], result, &
   y = [0.0_dp, 3.0_dp, 5.0_dp], exact = .true.)
if (result%status /= rank_ok) error stop "test failed"
print *, result%statistic, result%p_value
```

- `wilcox_exact(x,result,y,paired,alternative,mu,exact,conf_level)`:
  rank sum W or signed rank V; paired differences when requested; zero
  differences omitted; exact p-values include ties. Normal approximations
  follow exactRankTests and do **not** apply a continuity correction.
  `conf_level` requests exact shift intervals and a Hodges-Lehmann-style
  upstream estimate, including tied and shifted data.
- `ansari_exact(x,y,result,alternative,exact,conf_level)`: Ansari-Bradley
  scale test. Exact scale intervals/estimates require finite nonzero samples.
  The upstream asymptotic one-sided convention is preserved: its `less` and
  `greater` directions are reversed relative to its exact branch.
- `perm_test(x,result,y,paired,alternative,mu,exact,tolerance)`: raw-score
  permutation test; exact/noninteger default remapping follows upstream.
  Set `tolerance` explicitly to use the distribution lattice approximation
  instead of the test's default remapping.

Optional alternatives are the full strings `two.sided`, `less`, `greater`.
Tests default to two-sided. Omitted `exact` uses upstream size thresholds:
less than 50 for Wilcoxon/Ansari and at most 50 for permutation tests.
The result contains `statistic`, `p_value`, `point_probability`, `expected`,
`variance`, `z`, `estimate`, `conf_int(2)`, `has_conf_int`, `exact`, and `status`.
Rank-sum `statistic` and `expected` both use W's shifted scale. Point mass is
NaN for asymptotic tests; undefined standardized statistics/estimates are NaN.
`has_conf_int` records interval calculation, but endpoints can be NaN when
upstream's inversion cannot identify a confidence set.

Wilcoxon omits nonfinite values (paired omissions are joint). Ansari and
permutation tests omit NaNs but reject infinities. A wholly zero signed-rank
sample is invalid. Paired shape errors and invalid alternatives return a
status instead of stopping the program.

## Reusable distributions and simulation

```fortran
type(rank_distribution) :: law
integer :: status
call permutation_distribution([1.0_dp, 2.5_dp, 2.5_dp, 4.0_dp], 2, law, status)
if (status /= rank_ok) error stop "distribution failed"
print *, dperm(5.0_dp, law)
print *, pperm(5.0_dp, law, "greater")
print *, qperm(0.95_dp, law)
```

`permutation_distribution(scores,m,law,status,paired,factor,tolerance,max_cells)`
builds sorted `support` and `probability` arrays. Sampling is without
replacement. When `m == size(scores)`, the default is independent fair sign
flips of nonnegative scores; explicitly pass `paired=.false.` for a fixed-size
sample of the entire population. Otherwise `paired` cannot select sign flips.
Signed population scores are allowed for fixed-size sampling. Integer and
half-integer scores are exact lattices. Arbitrary finite real scores are rounded
with R's ties-to-even rule using a supplied `factor` or upstream's automatic
factor search; they are exact **for the rounded scores**, not generally for the
original reals. The default tolerance is 0.01. `factor` and `rounding_error`
(sum of absolute per-score errors) are reported in the distribution. An unmet
automatic error target is observable through that field rather than an R warning.

The shifted, greatest-common-divisor-compressed DP computes normalized
probabilities, avoiding overflow of binomial counts. Its default cap is
5,000,000 table cells (about 40 MB of DP storage), configurable with the optional
int64 `max_cells`. Scores whose scaled magnitude exceeds 2^50 are rejected
before integer conversion. A conservative bound rejects inputs with
size(scores)*max(abs(scores)) above the floating-point range, even if
cancellation could make a particular statistic finite. Oversized allocations are rejected rather than
silently switching to an approximation. Time is O(N*m*span) for fixed-size
sampling and O(N*span) for sign flips. Extremely small probabilities can
underflow to zero; such support points are omitted, so this is not an
arbitrary-precision probability implementation.

`dperm`, `pperm`, `qperm` are pure elemental functions over query values.
`pperm` defaults to `less`; both one-sided tails include equality. Unpaired
two-sided tails sum masses at least as far from the original-score expectation
as the observation. Paired two-sided tails double the smaller inclusive tail,
capped at one. `qperm` returns the smallest support point with CDF >= p;
p=0/1 returns the extreme retained support. Invalid probabilities/NaNs return
NaN. The raw arrays replace R's optional density data frame.

`rperm(scores,m,uniforms,draws,status)` draws fixed-size subset sums using a
matrix of independent [0,1) uniforms with shape `(m,B)`. It always samples
without replacement, including when m=N, matching the separate R `rperm`.
`simulate_distribution(scores,m,uniforms,law,status,paired)` instead builds
an empirical law, with the same sign-flip default as the exact builder.
Caller-managed uniforms make both routines pure and reproducible; neither
mutates a global RNG or reproduces R's exact random-number sequence.

## Scores and status

`cscores(y,scores,status,score_type,integer_scores,max_score)` supports Data,
Wilcoxon, NormalQuantile, AnsariBradley, Median, Savage, and ConSal.
`cscores_factor(levels,scores,status)` takes factor codes 1/2.
`cscores_surv(time,event,scores,status,integer_scores,max_score)` computes
log-rank scores with censoring and ties. `irank(x,ranks,status,order)` returns
maximum integer ranks and validates any supplied one-based sorting permutation.
Integer scaling preserves upstream's doubling shortcut for integer/half-integer
scores. Zero-denominator scaling returns a numerical status.

All numerical procedures are pure; scalar distribution queries are elemental.
Every maintained procedure uses the single `dp` imported from `r_kinds`.
Status constants: `rank_ok=0`, `rank_invalid=1`, `rank_limit=2`,
`rank_allocation=3`, `rank_numeric=4`. On failure, do not consume output arrays
or a test record's numerical fields. Only a law whose own `status` is `rank_ok`
can be queried. Some allocations for small auxiliary input-sized arrays still
use normal Fortran allocation failure handling.

## Deliberate boundaries

R formula evaluation, object construction, attributes, S3 dispatch, data sets,
messages/warnings, and partial-string argument matching are not reproduced.
Asymptotic confidence intervals and general permutation confidence sets remain
untranslated. Requesting a CI in asymptotic Wilcoxon/Ansari mode returns
`rank_invalid`. Ansari exact CI currently rejects zero observations and guards
ratio/evaluation overflows. Wilcoxon/Ansari CI grids are capped at one million
pairwise points; their direct step recomputation can be expensive. Undefined
upstream Wilcoxon estimates such as empty-extremum -Inf are normalized to NaN.

## Validation

`fpm test` runs 1005 deterministic checks: actual R 4.5.0/exactRankTests 0.8-37
fixtures for all test families and alternatives, tied/shifted exact intervals,
all seven score families in ordinary/integer mode, censored log-rank scores,
exact distributions/tails/quantiles, plus independent exhaustive subset/sign
enumeration, balanced supplied-uniform simulations, malformed-input checks,
and a 1,100-score case that would overflow unscaled combinatorial counts.
The checked-in fixture regenerates with `Rscript test/generate_reference.R`
from this directory; R is not needed to build or test the Fortran package.
Independent review additionally compared 250 enumerated distributions and
135 tied/shifted Wilcoxon cases against R.

## Translation coverage

Package status: **substantial**. **15 of 18 (83.33%)** exported computational R
functions have meaningful Fortran mappings. The fraction measures mapped
computational functions, not complete R compatibility. Generic entry points
and explicitly exported numerical S3 methods count separately, as required by
the translation convention. The three exported formula methods remain in the
denominator conservatively and are untranslated: `wilcox.exact.formula`,
`ansari.exact.formula`, `perm.test.formula`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `wilcox.exact` | `exact_rank_tests` | `wilcox_exact` | partial |
| `wilcox.exact.default` | `exact_rank_tests` | `wilcox_exact` | partial |
| `ansari.exact` | `exact_rank_tests` | `ansari_exact` | partial |
| `ansari.exact.default` | `exact_rank_tests` | `ansari_exact` | partial |
| `perm.test` | `exact_rank_tests` | `perm_test` | partial |
| `perm.test.default` | `exact_rank_tests` | `perm_test` | partial |
| `dperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `dperm` | substantial |
| `pperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `pperm`, `dperm` | substantial |
| `qperm` | `exact_rank_distribution` | `permutation_distribution`, `simulate_distribution`, `qperm` | substantial |
| `rperm` | `exact_rank_distribution` | `rperm` | substantial |
| `cscores` | `exact_rank_scores` | `cscores`, `cscores_factor`, `cscores_surv` | substantial |
| `cscores.default` | `exact_rank_scores` | `cscores` | substantial |
| `cscores.factor` | `exact_rank_scores` | `cscores_factor` | substantial |
| `cscores.Surv` | `exact_rank_scores` | `cscores_surv` | substantial |
| `irank` | `exact_rank_scores` | `irank` | substantial |

The table and counts match `fpm.toml`. Detailed boundaries appear in
[API_COVERAGE.md](API_COVERAGE.md).
