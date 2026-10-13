# Verification results

Validated 2026-10-11 with gfortran 14.2.0, FPM 0.13.0 and R 4.5.0 on Linux.
The live R comparison used Kendall 2.2.2 compiled from the pinned upstream
revision in NOTICE.md. The sibling rfortran-core came from target repository
revision cbd437d6d92f5481706fd8ea5a7829aaff7b6e24.

## Build, test and example runs

Passed:

- `fpm build`
- `fpm test`
- `fpm run --example kendall_demo`
- Debug build, tests and example with
  `--profile debug --flag '-Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace'`
- Release build, tests and example with `--profile release`
- `Rscript validation/run_reference.R`

No Windows execution is claimed. The source uses standard free-form Fortran,
portable FPM sibling dependencies, no platform-specific library calls and no
external BLAS/LAPACK links.

## Independent deterministic tests

- The two examples from the upstream Kendall manual: S = 34, the tie denominators,
  and the second example's exact variance 6709/33.
- Every permutation for n = 3 through 8 (46,230 permutations), independently
  counted by inversions. Every attainable absolute-score tail is tested against
  the exact method, including inclusive extremes and central probabilities.
- Three tie configurations, each enumerating all 7! index permutations:
  independent empirical mean/variance of S and direct non-tied pair counts
  verify all tie-variance terms and the denominator.
- Symmetry under variable exchange, sign reversal and two-sided p-value symmetry.
- Seasonal moment sums, unchanged missing time slots, constant seasons,
  partial final cycles and too-short complete seasons.
- Empty/short vectors, unequal lengths, all-missing data, invalid method/period,
  retained infinity, NaN-pair omission before finite validation and constants.
- Three-point approximate/exact inference, tied exact rejection, exact-size
  boundary and explicit rejection above 100, analytical 100-point exact tail,
  tiny normal tails, close binary64 values and finite endpoints near HUGE().

## R reference comparisons

The reproducible script checks:

- 240 association/Mann-Kendall cases against both the original package and a
  direct binary64 R reference, with ties and missing observations
- 31 seasonal cases, including unequal season lengths, shifted `ts` starts,
  missing observations and a constant seasonal component
- 96 exact-distribution cases: ordinary tails against R's `cor.test(...,
  method="kendall", exact=TRUE)` and monotonic extremes against `2/n!`

For the 120 ordinary paired cases, maximum error scaled by
`max(1, abs(upstream))` was `1.42862e-7`, consistent with the upstream binary32
implementation. Acceptance uses an absolute tolerance of 4e-7 plus a relative
tolerance of 2e-6 for upstream parity, and 3e-13 absolute/relative for the direct
binary64 reference. These tolerances do not apply to corrected degenerate,
three-point, close-value, or tiny-tail cases.

The analytical formula is intentionally used for extreme exact tails. In the
installed R 4.5.0, `cor.test(1:50, 1:50, method="kendall", exact=TRUE)` returned
approximately -2.13163e-14 due to tail cancellation; it is not a valid reference
for that extreme. The Fortran recurrence only sums nonnegative probabilities
and returns approximately 6.575898833e-65 there.

## Verified behavior differences

| Increasing sequence | Original Kendall p | Binary64 default p | Optional exact p |
| --- | ---: | ---: | ---: |
| n = 3 | 1, error 12 | 0.2962698714842864 | 0.3333333333333333 |
| n = 4 | 0.0894293785095215 | 0.0894293590289936 | 0.0833333333333333 |
| n = 5 | 0.0274863243103027 | 0.0274863361115104 | 0.0166666666666667 |
| n = 100 | 0 | 3.6371896876271862e-49 | 2.1430205762509339e-158 |

With `x = 1 + (1:5)*1e-10`, the upstream single-precision conversion collapses
all observations, returning D = 0 and a warning. This port preserves the order,
returns D = 10, and computes the ordinary five-point test.

Inspection and runtime agree: the upstream exact branch is never executed,
so its documentation's exact-p assertion is not used as a compatibility claim.
The default here retains the ordinary upstream normal formula; the exact method
is separately labeled an extension.

## Source and packaging checks

A conservative static audit passed: four maintained free-form Fortran files,
three distinct public R mappings, consistent counts/fraction/source paths,
explicit per-dummy INTENT and meaningful trailing FORD comments, one imported dp
kind, no disallowed real declarations/literals, no self-comparison NaN tests,
no semicolon-separated statements, no source lines over 132 characters,
no duplicate Fortran files and no copied dependency source.

Preserved upstream content is checksummed in upstream/SHA256SUMS. Original
fixed-form source is archival text only. FPM build products are temporary and
must be removed with `fpm clean --all` before packaging; this report does not
itself claim an archive was created or uploaded.
