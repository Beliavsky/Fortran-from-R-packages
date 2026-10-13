# PRROC: precision-recall and ROC curves in Fortran

A modern, free-form Fortran translation of the computational exports of
[CRAN PRROC 1.4](https://cran.r-project.org/package=PRROC), by Jan Grau and
Jens Keilwagen. PR AUC uses the package's **exact nonlinear piecewise
interpolation**, not trapezoidal integration of threshold precision values.
Weighted observations, binary and soft labels, exact score ties, the
Davis–Goadrich hard-label area, and optional reference curves are supported.

## Build and run

Place `PRROC/` beside the existing `rfortran-core/` directory in
[Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages).
From `PRROC/`:

```sh
fpm build
fpm test
fpm run --example
```

Only FPM and a modern gfortran are needed for these commands, including on
Windows. No R installation, external BLAS/LAPACK, system library linker flags,
or vendored dependency sources are required. The sibling path dependency
reuses `r_kinds:dp`, `r_ordering:r_order`, and the stable transforms in
`r_transforms`. This API was checked in repository revision
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.

## Example

```fortran
use prroc_api, only : dp, prroc_result, pr_curve, roc_curve, prroc_ok
type(prroc_result) :: pr, roc
integer :: status

call pr_curve([4.0_dp, 4.0_dp, 2.0_dp, 1.0_dp], pr, &
   weights_class0 = [0.2_dp, 0.9_dp, 0.7_dp, 0.1_dp], &
   curve = .true., rand_compute = .true., status = status)
if (status /= prroc_ok) error stop 'PR computation failed'
print *, pr%auc                 ! 0.559476125350028
print *, pr%random%auc          ! 0.475

call roc_curve([4.0_dp, 4.0_dp, 2.0_dp, 1.0_dp], roc, &
   weights_class0 = [0.2_dp, 0.9_dp, 0.7_dp, 0.1_dp], status = status)
if (status /= prroc_ok) error stop 'ROC computation failed'
print *, roc%auc                ! 0.650375939849624
```

The complete executable example is `example/classifier_curves.f90`.

## Public API

All computational procedures are `pure`; all floating-point values use
`real(dp)`, with `dp` re-exported from `r_kinds`.

- `pr_curve(scores_class0, result, scores_class1, weights_class0,
  weights_class1, sorted, curve, min_step_size, max_compute, min_compute,
  rand_compute, dg_compute, status, max_curve_points)`
- `roc_curve(scores_class0, result, scores_class1, weights_class0,
  weights_class1, sorted, curve, max_compute, min_compute, rand_compute,
  status, max_curve_points)`

Only `scores_class0` and `result` are required. The second argument is an
`intent(out)` `type(prroc_result)`; the remaining arguments are optional.
The `scores_class0` class is **positive/foreground**, following PRROC naming.
Higher scores predict the positive class. The following input forms work:

1. Separate positive and negative score vectors: provide both score arguments.
   Unspecified weights are ones when neither weight vector is present.
2. Binary labels or soft probabilities: supply all scores as `scores_class0`
   and positive labels/probabilities as `weights_class0`. Negative scores
   default to the same observations, and negative weights default to
   `1 - weights_class0`.
3. Independently weighted classes: provide both score and both weight vectors.
   Their lengths can differ; each weight vector must match its score vector.

Scores may contain infinities but not NaNs. Weights must be finite and
nonnegative. If the complementary default is used, it must be nonnegative
and its length must match the negative scores. `sorted=.true.` means both
score vectors are already in ascending order, and this is checked. Otherwise
stable sibling-library sorting retains the association between scores and
weights. R's lazy-default ordering is preserved: omitted second scores and
weights are resolved **after** sorting the first vector. In particular,
explicit second scores with an omitted second weight vector inherit the
complement of the already-sorted first weights. Supply both weight vectors
to avoid that unusual R convention for independently weighted classes.

### Results

`prroc_result` extends the public `prroc_summary` type. Both expose:

- `auc`: exact nonlinear PR AUC, or trapezoidal ROC AUC
- `auc_davis_goadrich`: the hard-label Davis–Goadrich PR AUC when eligible
  and requested; otherwise quiet NaN, including all ROC results
- `curve(:,:)`: an allocated real matrix, always with three columns; zero
  rows when not requested or unavailable

`prroc_result` additionally has allocatable `maximum`, `minimum`, and
`random` members, each of type `prroc_summary`. These are allocated only
when the corresponding `*_compute` flag is true. Do not access an
unallocated reference component.

`curve=.true.` returns `(recall, precision, threshold)` for PR, ordered in
**decreasing recall** like R; ROC returns `(false-positive rate,
sensitivity, threshold)` in increasing false-positive rate. Duplicate
endpoint and interpolated rows, and R's shifted ROC threshold-column
convention, are retained. The default is to return AUC only.

`min_step_size` is a positive finite increment in **original positive-weight
units**, not normalized recall units. Its R-compatible default is 1 for
unweighted input and `min(1, sum(weights_class0)/100)` when the first weights
are supplied. It changes PR sampling density, never the integrated AUC.
`dg_compute` defaults to true. Davis–Goadrich is eligible for distinct
unweighted class vectors or matched score vectors with complementary binary
weights; arbitrary integer or fractional weights do not enable it.

### Reference curves

`max_compute`, `min_compute`, and `rand_compute` default to false. Maximum
and minimum reproduce upstream's constructions: matched score vectors are
replaced by the first-class weights or their negatives; separate weighted
class vectors get class-indicator scores. These names denote the upstream
reference constructions; arbitrary independent, non-complementary weights
on matching score vectors do not guarantee a globally optimal ordering.
Random PR AUC is total positive mass divided by total mass. Random ROC AUC
is 0.5. Random references always contain two endpoint rows, even when
`curve` is false.

### Status and edge cases

`status` is optional but recommended. It reports:

| Named constant | Value | Meaning |
| --- | ---: | --- |
| `prroc_ok` | 0 | Requested numerical results computed |
| `prroc_invalid_input` | 1 | Empty overall input, invalid weights/shapes/step/limit, NaN score, or incorrect sorted promise |
| `prroc_degenerate` | 2 | No positive mass, or ROC also has no negative mass |
| `prroc_curve_limit` | 3 | Requested curve exceeds its row limit; AUCs are still available |

Errors 1 and 2 return quiet-NaN AUCs, an allocated zero-row curve, and no
references. Positive-only PR is well-defined here with AUC 1. ROC requires
positive mass in both classes. Zero-total-weight observations are removed
before threshold construction. Individual empty class vectors are allowed
when the remaining classes meet these requirements.

`max_curve_points` defaults to 1,000,000 rows per curve and must be at least
2. A curve exceeding the limit stays empty, while all requested AUCs and
any smaller reference curves remain available. Increase the limit or choose
a larger PR step if more rows are needed. This also bounds an otherwise
unbounded allocation for very large weights with the default step.

Separate class normalization and log-scale ratios preserve class imbalance
without overflow in total weights. Tiny within-class relative weights may
still underflow below the floating-point representable range. A separate
series protects the nonlinear integral from cancellation when a new segment
has extremely little mass. Do **not** use `-ffast-math`, `-Ofast`,
`-ffinite-math-only`, or equivalent flags: validation and unavailable-result
semantics depend on IEEE NaNs and infinities.

### Deliberate differences from R

- Typed results and status codes replace S3 lists, R errors/warnings,
  recycling, and presentation methods.
- Nonfinite/negative weights, shapes and sorted promises are validated for
  both functions. PRROC's ROC implementation does not consistently validate
  weights.
- Zero-mass thresholds are removed, preventing artificial undefined initial
  precision values and extra curve rows. No-positive data return a clear
  degeneracy status rather than propagated NaNs/warnings.
- The random PR denominator uses `sum(w0) + sum(w1)`. Upstream's
  `sum(w0 + w1)` recycles unequal-length weight vectors and can give an
  incorrect baseline. This correction is independently tested.
- Unused all-NA padding in R's random curve matrices is omitted.
- Stable scaling/logarithms avoid overflow and improve small-increment
  accuracy, so extreme cases can remain finite where R loses precision.
- Curve allocation is explicitly bounded. Plotting, `ggprroc`, printing,
  colors, and R object construction are omitted.

## Translation coverage

Package status: **substantial**. **2 of 2 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R interface or edge-case compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `pr.curve` | `prroc_api` | `pr_curve` | substantial |
| `roc.curve` | `prroc_api` | `roc_curve` | substantial |

`ggprroc` is presentation-only and excluded from the denominator. Registered
S3 plotting/printing methods and internal helpers are also excluded. No
computational export remains unmapped. See `API_COVERAGE.md` and the matching
`extra.translation` section in `fpm.toml`.

## Validation

The checked-in generated test was produced by actually executing the
preserved PRROC 1.4 source under R 4.5.0. It contains **762 R-reference
assertions across 43 deterministic cases**, including complete observed and
reference curve matrices, soft labels, ties, sorted inputs, lazy defaults,
custom/default curve steps and seeded randomized cases. Comparisons use
absolute tolerances of 2e-12 for scalar areas and 3e-12 for curves. Comparisons
explicitly omit changed zero-weight curve rows and the corrected unequal-
length random denominator; they are not silently substituted R outputs.

The independent suite adds **85 assertions** for weighted pairwise ROC,
score-reversal symmetry, exact nonlinear PR identities, scale invariance at
1e250 and 1e-250, 500-order class imbalance, a 5e-21 area, overflowing raw mass
totals, zero weights, absent classes, input validation, and allocation limits.

To regenerate the R fixture, run from this directory:

```sh
Rscript tools/generate_reference.R
```

The preserved source is loaded directly; no installed PRROC or `rlang` is
needed for regeneration. Normal FPM builds/tests do not execute R.

Validation uses gfortran 14.2.0 and FPM 0.13.0 in default, checked debug
(`-g -O0 -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Wall -Wextra`),
and release profiles. Each profile builds, tests and runs the example.
Clean between profiles to avoid stale executable cache behavior:

```sh
fpm clean --all
fpm build --profile release
fpm test --profile release
fpm run --example --profile release
fpm clean --all
```

## License and attribution

The translation and retained upstream reference material are **GPL-3.0-only**,
consistent with PRROC's `License: GPL-3`; see `LICENSE`, `NOTICE.md`, and
`upstream/PROVENANCE.json`. The independent sibling `rfortran-core` dependency
has its own MIT license and is not copied here. Preserve all these notices
when redistributing. Scientific references and the original citation metadata
are in `NOTICE.md` and `upstream/CITATION`.
