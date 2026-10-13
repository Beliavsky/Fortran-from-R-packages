# outliers

Modern free-form Fortran/FPM translation of CRAN **outliers 0.15**, by Lukasz Komsta.
Implements Grubbs, Dixon, chi-square and Cochran tests; their critical values and
probabilities; endpoint selection/removal; five outlier-score families; and the
upstream table-interpolation routine. GPL-2.0-or-later; see LICENSE and NOTICE.md.

## Build and run

Place this directory as `outliers/` at the root of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
next to the existing `rfortran-core/` translation. Then run from `outliers/`:

```sh
fpm build
fpm test
fpm run --example demo
```

Only gfortran and FPM are needed, including on Windows. There are no system
BLAS/LAPACK links, R runtime calls, network downloads, or vendored dependencies.
The sibling FPM dependency supplies normal/t/F/chi-square distributions, sorting,
quantiles, MAD, descriptive statistics and the single real kind `r_kinds::dp`.
The tested repository revision is recorded in NOTICE.md. All maintained routines
are pure, with elemental distribution functions accepting conformable arrays.
Do not use fast-math, Ofast, finite-math-only, or other flags that discard IEEE
NaN/infinity semantics. Default and release builds are supported.

## Quick use

```fortran
use outliers_api, only : dp, grubbs_test, outlier_test_result
real(dp) :: x(6)
type(outlier_test_result) :: result
x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 20.0_dp]
result = grubbs_test(x, type = 10, two_sided = .true.)
if (result%status == 0) print *, result%statistic, result%p_value
```

The full API is available through `outliers_api`, including `dp`. Numerical
distribution functions are also available directly from `outliers_distributions`.
Every argument has a trailing FORD comment in the source.

## API and result shapes

* `outlier(x, opposite)` returns one endpoint for a vector or one per matrix
  column. Selection compares distances from the sample mean; exact ties select
  the maximum. `outlier_mask(x, opposite)` selects every equal endpoint.
* `rm_outlier(x, fill, median, opposite)` removes every selected tie, preserving
  order, or replaces them with the mean/median of the remaining sample.
  Vector output is allocatable. Matrix output is an array of `outlier_column`,
  each with allocatable `%values`, so unequal post-removal lengths are supported.
* `scores(x, type, prob, lim)` preserves vector/matrix shape and scores each
  column separately. Types are exact names `z`, `t`, `chisq`, `iqr`, `mad`.
  Absent/NaN `prob` returns scores; `prob=1` returns CDF values; other probabilities
  return real 0/1 threshold flags, preserving NaN rather than forcing it into a
  Fortran logical. `prob=0` implements the source maximum-z rule except for
  chi-square scores. IQR uses `lim`, ignoring `prob`. The MAD constant is 1.4826.
* `chisq_out_test(x, variance, opposite)`, `dixon_test(x, type, opposite, two_sided)`
  and `grubbs_test(x, type, opposite, two_sided)` return `outlier_test_result`.
  `%statistic` is chi-square, Q or G, except type 20 Grubbs where it is U.
  `%u` is the residual/total sum-of-squares ratio for every Grubbs variant.
  `%selected` contains one or two tested values. `%n` counts cleaned observations;
  `%variant`, `%lower`, `%p_value`, and `%status` describe the computation.
  Status 0 means success; status 1 means invalid data/domain, with NaN numerical
  results (a statistic may remain available if only its probability is invalid).
  Unused numerical result fields are NaN, and unused `%groups` is empty.
* `cochran_test(variances, counts, inlying)` accepts group variances and real-valued
  group counts. `cochran_test(matrix, inlying)` takes equal-size groups by column.
  `cochran_groups(x, groups, inlying)` accepts unequal-size groups with consecutive
  integer labels 1 through k. `%statistic` is C, `%df` the mean group size,
  `%n` the number of groups, `%group` the first selected group and `%groups` all ties.
  Every group must contain at least two observations.
* `qdixon(p,n,type,rev)`/`pdixon(q,n,type)` and
  `qgrubbs(p,n,type,rev)`/`pgrubbs(q,n,type)` accept integer sample size.
  `qcochran(p,n,k)`/`pcochran(q,n,k)` accept real mean group size n>1 and integer k>=2.
  Quantile/CDF procedures are elemental. Dixon probabilities are upper-tail
  values, unlike Grubbs/Cochran CDFs. Follow the R conventions, not generic naming
  assumptions. Grubbs type 20 uses lower-tail U probabilities.
* `qtable(p,probs,quants)` accepts scalar/vector evaluation coordinates and paired
  coordinate/value vectors, sorted internally. It uses the same four-node choice
  as R for cubic interpolation and two endpoint nodes for linear extrapolation.
  A duplicate endpoint pair produces its intercept-only mean; duplicate nodes
  inside the chosen cubic stencil return NaN, matching an R fitting error.

Defaults follow R: Dixon type 0 and two-sided true; Grubbs type 10 and two-sided
false; other logical flags false. Tests use R's doubled-and-reflected two-sided
p-value transformation. No formula parsing, data-frame/S3 construction, attributes,
warning strings, partial argument-name matching, or R vector recycling is provided.
Cochran formulas are replaced by typed grouped observations.

## Domain and fidelity notes

* NaNs represent missing R NA/NaN. Endpoint selection and the three individual
  observation tests omit them. Logical endpoint masks use false at missing
  positions, whereas R can return NA there. Removal preserves missing positions;
  replacement uses R's non-NA-removing remaining-data mean/median. An all-missing
  endpoint is NaN; removing from an all-missing vector leaves it unchanged.
* Scores propagate missing data. All data-based routines reject nonfinite
  nonmissing input; some R IQR/MAD edge cases with infinities therefore differ.
  Empty/degenerate tests return status 1 instead of raising R errors or producing
  accidental finite p-values. Grubbs requires n>=3 and positive finite variance;
  constant input is invalid rather than reproducing upstream's spurious p=0.
  A three-point type 11 Grubbs test can return G/p with undefined U (NaN).
  Chi-square tests require positive variance; a singleton is supported only
  when an explicit positive reference variance is supplied.
* Dixon type 10 supports n=3..30; types 11/20 n=4..30; types 12/21 n=5..30;
  type 22 n=6..30. Automatic type 0 chooses 10 through n=7, 11 through n=10,
  21 through n=13, then 22. No larger-sample extrapolation is invented.
* Grubbs type 20 has usable table data only for n=4..30. Types 10/11 use shared
  Student t distributions, with forward p in [0,1); p=1 retains R's undefined
  endpoint as NaN. For type 11 reverse mode, the source's clipping at forward
  probabilities 2e-16 and .9999 and its .001 lower root bracket are retained.
  The narrow bracket gap returns NaN rather than an R uniroot error. Within
  the valid bracket, an algebraic inversion replaces the approximate root solve.
  Comparisons within 128 machine epsilons times the statistic scale recognize
  boundary equality, avoiding backend-rounding changes to clipping. At the lower
  discontinuity, R can produce zero or a bracket error one ulp apart; equality
  returns NaN here. The probability thresholds themselves are unchanged.
* Finite coordinates outside tabulated probabilities use the original endpoint
  extrapolation, then Dixon/Grubbs-table outputs are clamped to [0,1]. The source
  tables are empirical approximations, not exact mathematical distributions.
  Forward and reverse cubic interpolations are not exact inverses between nodes.
* **Upstream anomalies are deliberately preserved.** Dixon type 21, n=7, p=.7
  contains 350; type 22, n=6, p=.4 contains 830. Other nonmonotone entries and
  the duplicated .009 in type 12, n=29 also remain unchanged. This can produce
  clipped, nonmonotone or undefined results near anomalous entries. Use scientific
  judgment; the translation does not silently correct the reference tables.

## Validation

`test/*.dat` are deterministic numerical fixtures generated by actual R and
outliers 0.15 using `tools/generate_references.R`. They cover every Dixon sample
size/type, all table nodes, interpolation/extrapolation, automatic types, both
directions, all Grubbs variants, extreme probabilities, Cochran parameters, test
options, missing values, score modes and small samples. R need not be installed to
run `fpm test`. Additional independent checks cover cubic polynomials, inverse
relations, exact published table values, anomalies, ties, matrices, ragged removal,
replacement and Cochran group shapes. See VALIDATION.md for executed checks.

To regenerate tables: `python tools/generate_tables.py`. To regenerate fixtures
with installed R outliers 0.15: `Rscript --vanilla tools/generate_references.R .`.
Run only one generator at a time. Preserved upstream sources are reference material
only and are not dependency implementations.

## Translation coverage

Package status: **substantial**. **14 of 14 (100%)** exported computational R
functions have meaningful mappings. This measures mapped computational functions,
not complete R compatibility. Every function exported by `exportPattern(".*")`
is included, including the computational helper `qtable`. The exact same mappings
and status values are recorded in `fpm.toml`; API_COVERAGE.md adds per-function notes.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `chisq.out.test` | `outliers_api` | `chisq_out_test` | substantial |
| `cochran.test` | `outliers_api` | `cochran_test`, `cochran_groups` | substantial |
| `dixon.test` | `outliers_api` | `dixon_test` | substantial |
| `grubbs.test` | `outliers_api` | `grubbs_test` | substantial |
| `outlier` | `outliers_api` | `outlier`, `outlier_mask` | substantial |
| `pcochran` | `outliers_distributions` | `pcochran` | complete |
| `pdixon` | `outliers_distributions` | `pdixon` | substantial |
| `pgrubbs` | `outliers_distributions` | `pgrubbs` | substantial |
| `qcochran` | `outliers_distributions` | `qcochran` | complete |
| `qdixon` | `outliers_distributions` | `qdixon` | substantial |
| `qgrubbs` | `outliers_distributions` | `qgrubbs` | substantial |
| `qtable` | `outliers_distributions` | `qtable` | complete |
| `rm.outlier` | `outliers_api` | `rm_outlier` | substantial |
| `scores` | `outliers_api` | `scores` | substantial |
