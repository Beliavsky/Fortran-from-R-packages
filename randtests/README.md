# randtests

Modern free-form Fortran translation of **randtests 1.0.2**, by Frederico Caeiro
and Ayana Mateus. Six nonparametric randomness tests, runs distributions, exact
Bartels distributions, and ordered permutations are provided.

## Build and use

Place this directory at the repository root beside `rfortran-core`:

```text
Fortran-from-R-packages/
  rfortran-core/
  randtests/
```

From `randtests/`:

```sh
fpm build
fpm test
fpm run --example randomness_demo
```

Only FPM, a Fortran 2008 compiler, and the sibling `rfortran-core` dependency are
required. There are no system BLAS/LAPACK links, network dependencies, copied
shared modules, or R runtime dependencies. The manifest is portable to Windows
with gfortran and FPM; validation here used Linux, gfortran 14.2.0 and FPM 0.13.0.
Do not use `-ffast-math`, `-Ofast`, or `-ffinite-math-only`: IEEE NaN/Inf semantics
are part of this API.

```fortran
use randtests, only: dp, randomness_result, runs_test
real(dp) :: x(6) = [1.0_dp, 4.0_dp, 2.0_dp, 6.0_dp, 3.0_dp, 5.0_dp]
type(randomness_result) :: result
result = runs_test(x, pvalue='exact')
if (result%status /= 0) error stop 'Runs test failed'
print *, result%statistic, result%p_value
```

`use randtests` re-exports the complete public API and the shared `dp` kind.
All deterministic procedures are `pure`; scalar distribution procedures are
`pure elemental`. `rruns` alone changes the runtime RNG state.

## API and numerical contract

- Test functions take `real(dp) :: x(:)`, optional `alternative`, and return
  `randomness_result`. Alternatives are `two.sided` (default), `left.sided`,
  `right.sided`, or aliases `t`, `l`, `r`. All tests consistently omit NaNs and
  reject infinities. NaN omission in Cox–Stuart and turning-point tests extends
  the inconsistent upstream behavior.
- Result fields: `statistic`, `p_value`, `observed`, `mean`, `variance`, `rvn`,
  `n`, `n1`, `n2`, `status`, and `message`. `observed` is the unnormalized count
  (or NM for Bartels); `rvn` is populated for Bartels only; `n1`/`n2` are
  populated for runs only. Unused successful fields are zero. `n` is the
  post-preprocessing sample size, except Cox–Stuart uses the number of unequal
  pairs. `mean`/`variance` refer to RVN for Bartels and to `observed` otherwise.
- Status constants: `rt_success=0`, `rt_invalid=1`, `rt_insufficient=2`,
  `rt_degenerate=3`, `rt_limit=4`. Error results have NaN real fields and a
  diagnostic message; no routine prints, stops, or raises an R-style exception.
  Always check `status` before reading test results.
- `runs_test` accepts optional finite `threshold` (default: median after NaN
  removal), and `pvalue='normal'` or `'exact'`. Threshold-equal observations are
  removed. Both classes must remain. When n1=n2=1, exact mode returns p=1 and
  statistic=0 by convention; normal mode reports degeneracy.
- `bartels_rank_test` accepts `pvalue='normal'`, `'beta'`, `'exact'`, or `'auto'`.
  Auto chooses exact for n<=9, beta for n<=100, otherwise normal. Exact mode
  rejects ties and n>9; auto consequently rejects tied samples of size <=9.
  Supply `'beta'` or `'normal'` explicitly for small tied samples. Average ranks
  are used with ties, but the original untied approximation variance is retained.
  Constant sequences are rejected. These normal/beta approximations are not a
  tie-adjusted exact test. `rank_test` likewise retains upstream's unadjusted
  pair-count statistic and variance when ties occur.
- `druns(x,n1,n2[,log_probability])`, `pruns(q,n1,n2[,lower_tail,log_probability])`,
  and `qruns(p,n1,n2[,lower_tail,log_probability])` are scalar/elemental real
  functions. Counts must be positive integers and total <=100000. Invalid
  parameters produce NaN, zero log probabilities produce -Inf. Noninteger PMF
  points have zero mass; CDF cutoffs are floored. Infinite CDF cutoffs are valid.
  The upper tail is strictly P(R>q). The quantile is the smallest support value
  whose lower CDF is at least the requested lower probability. A relative
  log-probability tolerance of 64 machine epsilons protects discrete boundaries
  against rounding in exp/log conversions, with boundary quantiles equal to the minimum/maximum support.
- `rruns(values,n1,n2[,status])` fills a caller-sized integer array using
  `random_number`; invalid counts yield -1 elements and status=1. The caller
  controls `random_seed`. Repeated seeds reproduce results on the same runtime;
  neither cross-compiler sequences nor parity with R's RNG are promised.
- `dbartelsrank(x,n[,log_probability])` and
  `pbartelsrank(q,n[,lower_tail,log_probability])` enumerate untied ranks 1:n for
  1<=n<=9. `bartels_histogram(n,histogram,status)` exposes exact integer counts
  with lower bound zero. Enumeration uses a small histogram and recursion, not
  a factorial-sized matrix, but costs O(n!) time and is recomputed per call.
- `permut(x,matrix,status[,m])` returns one ordered selection per row. Duplicate
  values remain distinct source positions. Default m=size(x); 0<=m<=size(x).
  m=0 returns a 1-by-0 matrix. At most 1000000 output elements are permitted.
  Invalid input gives status=1, resource limits status=2, and a 0-by-0 matrix.
  Row order follows input-position lexicographic traversal and may differ from
  upstream's hard-coded small-n ordering.
- `permut_apply(x,fun,values,status[,m])` applies a pure scalar-real callback to
  every row. Declare `procedure(permutation_callback), pointer :: fun`, associate
  it with a compatible pure function, and pass the pointer. Null pointers give
  status=1. Callback scalar results (including any NaN they produce) pass through.
  A pointer is used so every dummy, including this procedure dummy, has explicit
  `intent(in)`. R variadic arguments and vector-valued callbacks are omitted.
- `randtests_aux(x,method)` implements `bartels` (sum of squared adjacent
  differences) and `difference.sign` (number of increases). Empty/singleton
  inputs yield zero. Invalid method, nonfinite data, or an unrepresentable
  squared sum gives NaN.

R lists, S3 `htest` construction, names, warnings, automatic recycling, language
objects, datasets, printing, and plotting are not translated. Elemental calls
accept scalar arguments or conformable arrays, following Fortran rules.

## Deliberate corrections to upstream 1.0.2

This package implements the documented statistical definitions rather than
reproducing known upstream defects:

1. Bartels exact tests use NM=sum(diff(rank(x))**2). The R exact branch passes
   sum(rank(x)**2) instead. Exact right-tail tests here use P(NM>=observed), and
   two-sided tests use min(1,2*min(P(NM<=observed),P(NM>=observed))). The R branch
   incorrectly uses its lower CDF at observed-1 as the right tail.
2. `pruns(log.p=TRUE)` upstream cumulatively sums log masses; here the logarithm
   of the summed probability is returned, using log-sum-exp. Upper tails are
   summed directly. Real CDF cutoffs use floor, consistent with P(R<=q), whereas
   upstream rounds to nearest integer.
3. `qruns` honors `lower_tail` and `log_probability`, unlike upstream, and gives
   defined errors for invalid probabilities. Endpoints are handled explicitly.
4. Permutation callbacks run for every selection, including singleton and
   two-element cases skipped by upstream. This also corrects exact Bartels
   distribution functions for n=1 and n=2.
5. Two-sided discrete p-values are capped at one. The upstream Cox–Stuart and
   runs implementations can return values above one; all-tied Cox pairs give
   p=1 here. NaN handling, invalid-method diagnostics, finite-input checks, and
   bounded enumeration are also intentionally stricter.

For example, ranks [1,2,3] have NM=2, exact Bartels left p=1/3, right p=1, and
two-sided p=2/3. For n1=2,n2=3, P(R<=3)=1/2 and its log is log(1/2).
These corrected branches are tested against independent enumeration, not against
faulty R outputs. See `API_COVERAGE.md` and `validation/README.md`.

## Translation coverage

Status: **substantial**. **14 of 14 (100%)** exported computational R functions
have a meaningful mapping. This fraction measures mapped computations, not full
R compatibility. The conservative denominator includes `randtests.aux`, although
its Rd page is labeled internal: `NAMESPACE` uses `exportPattern(".")`, and the
function is exported and documented. `permut` is partial because of callback/type
restrictions; all other mappings are substantial for the differences above.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bartels.rank.test` | `randtests` | `bartels_rank_test` | substantial |
| `cox.stuart.test` | `randtests` | `cox_stuart_test` | substantial |
| `difference.sign.test` | `randtests` | `difference_sign_test` | substantial |
| `rank.test` | `randtests` | `rank_test` | substantial |
| `runs.test` | `randtests` | `runs_test` | substantial |
| `turning.point.test` | `randtests` | `turning_point_test` | substantial |
| `druns` | `randtests_distributions` | `druns` | substantial |
| `pruns` | `randtests_distributions` | `pruns` | substantial |
| `qruns` | `randtests_distributions` | `qruns` | substantial |
| `rruns` | `randtests_distributions` | `rruns` | substantial |
| `dbartelsrank` | `randtests_distributions` | `dbartelsrank` | substantial |
| `pbartelsrank` | `randtests_distributions` | `pbartelsrank` | substantial |
| `permut` | `randtests_permutations` | `permut`, `permut_apply` | partial |
| `randtests.aux` | `randtests_permutations` | `randtests_aux` | substantial |

## Provenance and license

GPL-2.0-or-later, preserving upstream's GPL (>=2) grant. See `LICENSE`,
`LICENSE-GPL-3`, and `NOTICE.md`. The exact CRAN mirror source snapshot is retained
under `upstream/` for review and reproducible reference generation. The original
`permut` credits `gtools::permutations`; that attribution is preserved. Shared
MIT-licensed numerical helpers are used through the sibling dependency, never
copied into this package.
