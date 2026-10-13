# monobin

Modern free-form Fortran translation of **monobin 0.2.4**, Andrija Djurovic's
monotonic binning algorithms for credit rating models. Binary and continuous
responses are supported by all seven algorithms. GPL-3.0-or-later.

## Build and run

Place this directory at the repository root beside **rfortran-core** and
**Hmisc**, then run:

```sh
fpm build
fpm test
fpm run --example binning
```

No R installation, BLAS, LAPACK, or external numerical libraries are needed for
Fortran builds, including gfortran/FPM builds on Windows. Dependencies are sibling
FPM paths and are deliberately not included in this package. Do not compile with
`-ffast-math`, `-Ofast`, `-ffinite-math-only`, or similar flags: NaN and infinity
semantics are part of the API. `dp` is imported from `r_kinds` and re-exported.
Native Windows execution remains a platform check separate from Linux testing.

## API

```fortran
use monobin, only : dp, bin_options, bin_result, iso_bin, apply_bins, bin_ok
type(bin_options) :: options
type(bin_result) :: fitted
integer, allocatable :: groups(:)
real(dp), allocatable :: estimated_risk(:)
! x and y are same-length real(dp) arrays.
options%direction = 1
call iso_bin(x, y, fitted, options)
if (fitted%status /= bin_ok) error stop 'Invalid or degenerate data'
call apply_bins(fitted, new_x, groups, estimated_risk)
```

Each algorithm has the interface
`call NAME_bin(x, y, result, options, special)`, with optional `options` and
`special`. `fit_bins(x,y,method,result,options,special)` dispatches on the exact
three-letter names `pct`, `iso`, `woe`, `cum`, `sts`, `ndr`, and `mdt`.
All maintained computational procedures are pure; eligible scalar helpers are elemental.

### Options

| Field | Default | Meaning |
|---|---:|---|
| `g` | 0 | Automatic default: 15 for pct/cum and 50 for mdt; positive values override it |
| `direction` | 0 | Automatic Spearman sign; 1 increasing or -1 decreasing |
| `target` | 0 | Auto-detect 0/1; 1 requires binary, 2 forces continuous |
| `separate_special` | false | Pool all special values, or keep distinct special categories |
| `woe_trend` | true | pct continuous trend check uses WoE; false checks response means |
| `min_pct_obs` | .05 | iso/woe/sts/ndr/mdt minimum fraction, with the original 30-row floor |
| `min_avg_rate` | .01 | Minimum binary fraction of all events, rounded up with floor 1; continuous minimum mean |
| `woe_gap` | .1 | woe requires strictly larger adjacent evidence gaps |
| `p_value` | .05 | sts/ndr retain only contrasts with strictly smaller p values |

Fields unused by an algorithm have no algorithmic effect, but option-domain
validation is shared. Fractions for observations must be in [0,1], minimum rate
and WoE gap nonnegative, and significance in (0,1]. Zero `g` means default and
negative `g` is invalid. No recycling is performed.

### Results and prediction

- `status`: `bin_ok=0`, `bin_invalid=1`, `bin_no_data=2`,
  `bin_constant_x=3`, `bin_constant_y=4`. The two constant-case diagnostics
  apply to complete cases, as upstream. On failure, `bins` is empty and other
  output must not be treated as a fit.
- `bins(:)`: special bins first, followed by nonoverlapping complete bins in
  ascending predictor order. Each `bin_summary` includes count, `y_sum`,
  `y_mean`, `x_mean`, `x_min`, `x_max`, `woe`, `iv`, distributions, and response
  sufficient statistics. Centered `y_m2` is used for stable inference.
- `membership(:)`: one-based indices into `bins`, with the **original input
  length**. Missing-target positions are 0. `retained(:)` gives the original
  indices of retained observations; `membership(retained)` corresponds to R's
  shorter `x.trans` vector.
- `breaks(:)`: lower observed predictor boundaries for complete bins 2 through K.
  `direction` records the supplied or raw-Spearman direction. With continuous
  pct/WoE checks and negative response totals, WoE direction need not equal
  response-mean direction.
- `dist_good/dist_bad`: binary good/bad proportions. For continuous responses,
  these fields contain response-sum share and observation share respectively.
  WoE is log(first share / second share); IV is their difference times WoE.
- `apply_bins(model,x,membership,response,woe)` returns optional fitted response
  means and WoE. Complete intervals are left-inclusive, with unbounded tails;
  new points in gaps between observed predictors are assigned by `breaks`.
  Unseen separate special categories return index 0 and NaN predictions. With
  pooled special cases, any declared/nonfinite special value maps to the fitted
  pooled special bin if it exists.

NaN response values are discarded. Infinite responses are rejected. Nonfinite
predictors are always special cases; extra finite sentinel values may be passed
in `special(:)`. In separate mode, categories appear in first-observation order,
rather than R's formatted-string ordering. R's NA and NaN categories are one IEEE
NaN category in Fortran. Special predictor summaries propagate NaNs.

## Numerical behavior and source compatibility

The translation implements the actual seven methods, not seven aliases of a
single isotonic smoother:

1. **pct** decreases the requested Hmisc-style equal-count group count until
   the required strict mean/evidence trend holds.
2. **iso** performs weighted PAVA after combining equal predictors. R's
   `isoreg` sorts tied predictors by descending signed response; this yields the
   same atomic tied-predictor pools. Count/rate corrections scan in increasing
   response order, including reversed predictor order for decreasing fits.
3. **woe** repeatedly merges the smallest adjacent evidence gap after iso.
4. **cum** repeatedly selects the first maximum cumulative response mean in
   the required direction. Equal-rate adjacent bins can remain distinct.
5. **sts** uses one-sided pooled two-proportion or Welch t-test comparisons.
6. **ndr** backward-deletes the weakest nested-dummy contrast. Saturated
   nested models are computed in the equivalent independent group-intercept
   parameterization: logistic IRLS/Wald tests for binary targets, common residual
   variance and t tests for continuous targets. Logistic fitting uses R's
   25-iteration and 1e-8 deviance convergence defaults. Explicit final cuts use
   the shared `Hmisc::cut2_explicit` API.
7. **mdt** preserves upstream breadth-first threshold updates, the original
   restrictive candidate-index limits, and extraction of cutpoints attached to
   final leaves. Gini and SSE selection share the equivalent squared-error
   ordering for binary responses.

Hmisc's automatic `cut2(g)` uses upper ECDF endpoints, consumes entire tied runs,
and suppresses repeated cuts. It is **not** ordinary interpolated quantiles.
The existing Hmisc translation exposes only explicit cutpoints, so the small
missing automatic-group adapter is implemented here with attribution. Numerical
sorting, average ranks and distribution CDFs reuse rfortran-core. No translated
dependency package source is vendored.

### Deliberate differences and upstream edge behavior

- Formatted labels, data-frame/S3 construction, R warning text and global-option
  mutation are omitted. Numeric group IDs are authoritative, avoiding rounded
  label collisions. There is no exported upstream prediction method; the typed
  apply operation makes the label-implied intervals reusable and fills gaps and
  out-of-range values that `slice.variable` leaves unchanged.
- No WoE pseudocount is added. Zero binary cells retain infinite evidence and IV;
  continuous nonpositive log arguments retain NaN. pct's continuous default can
  raise `missing value where TRUE/FALSE needed` upstream when a candidate has
  NaN WoE. Here such a candidate fails the trend check and grouping continues
  to one bin. The case is explicitly covered by a failure/extension test.
- `cum.bin.aux` has an early return without its `type` column when its reversed
  prebins are strictly increasing, including some contradictory forced-trend
  inputs. Here cumulative pooling still runs and returns a consistent typed fit.
- NDR's selected cutpoint at the observed maximum is an **upper endpoint** in
  upstream `cut2(cuts=...)`, so it merges the final two model groups. This
  surprising behavior is preserved through the shared explicit-cut routine and
  has regression tests.
- Finite precision can change decisions exactly at equal objective values or
  a requested p-value threshold. Inference uses centered sums of squares rather
  than cancellation-prone raw second moments. Near-singular/constant Welch
  contrasts are handled deterministically (equal means p=.5 one-sided, distinct
  zero-variance means p=0), whereas R may warn, return NaN, or stop.
- Special NA/NaN distinctions and character-label ordering cannot be represented
  identically by a numeric IEEE-only interface. Other counts, cutpoints, fitted
  memberships, means, evidence and information values are parity-tested.

## Tests

`fpm test` runs **125 successful genuine upstream-R fixtures**, generated by the actually installed monobin 0.2.4 exported namespace under
R 4.5.0 with dplyr 1.2.1 and Hmisc 5.3-0 `cut2`. The installed-namespace
fixtures are byte-identical to the independent original-source execution.
It checks exact bin counts, memberships and special counts, numeric boundaries,
response/predictor summaries, WoE, IV and STS adjacent p-values. Numeric tolerance is 3e-9 relative
with a unit absolute scale. Fixtures span binary/continuous responses, both
forced directions, ties, irregular predictors, 1,000 requested groups, rare
events, pure cells, missing targets, pooled/separate special values, negative
responses and nondefault minimums/gaps/significance.

Separate independent tests check assignment and prediction consistency,
response/count conservation, calibration, normalized distributions, ordered
nonoverlapping bins, forced monotonicity, boundaries and tails, unseen special
categories, and invalid/degenerate inputs. `test/regenerate.R` regenerates the
fixtures using a genuine installed monobin 0.2.4; R is needed only to regenerate.
One additional mixed-sign pct/WoE case records the upstream error and is covered
by the documented Fortran continuation behavior rather than a fabricated golden.

## Translation coverage

Package status: **substantial**. **7 of 7 (100%)** exported computational R
functions have meaningful translations. This fraction measures mapped
computational functions, not complete R interface compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `pct.bin` | `monobin` | `pct_bin` | substantial |
| `iso.bin` | `monobin` | `iso_bin` | substantial |
| `woe.bin` | `monobin` | `woe_bin` | substantial |
| `cum.bin` | `monobin` | `cum_bin` | substantial |
| `sts.bin` | `monobin` | `sts_bin` | substantial |
| `ndr.bin` | `monobin` | `ndr_bin` | substantial |
| `mdt.bin` | `monobin` | `mdt_bin` | substantial |

Internal aggregation, fitting, evidence and assignment helpers are implemented
as necessary but not double-counted. The `gcd` dataset, examples and presentation
helpers are excluded from the computational denominator. See `API-COVERAGE.md`,
`NOTICE.md`, `LICENSE`, and the matching manifest entries.
