# trend

Modern free-form Fortran/FPM translation of Thorsten Pohlert's R package
**trend 1.1.8**: trend tests, Sen slopes, change-point detection, and randomness
tests. All numerical public procedures are **pure**. No R installation is
needed to build, test, or run the Fortran package.

## Build

Extract/place this directory as `trend/` at the repository root, alongside
`rfortran-core/`, `Kendall/`, `randtests/`, and `exactRankTests/`:

```text
fpm build
fpm test
fpm run --example trend_demo
```

Dependencies are sibling paths in `fpm.toml`; they are not bundled or vendored.
The dependency closure uses no external BLAS/LAPACK, system linker flags, C
runtime extension, or R runtime. It is intended for FPM and gfortran on Windows
as well as other platforms. Linux gfortran 14.2/FPM 0.13 was actually tested;
a native Windows run has not been performed.

Never enable `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent flags:
NaN, infinity, and missing-value policies depend on IEEE semantics.

## Quick example

```fortran
use trend_module, only : dp, trend_result, trend_ok, mk_test, sens_slope
real(dp) :: x(5) = [1.0_dp, 2.0_dp, 1.5_dp, 3.0_dp, 4.0_dp]
type(trend_result) :: fit
call mk_test(x, fit, alternative = 'greater')
if (fit%status == trend_ok) print *, fit%statistic, fit%p_value
call sens_slope(x, fit)
print *, fit%estimate, fit%conf_int
```

## API

Import `trend_module`; it re-exports the shared `dp` kind. All arrays use
`real(dp)`. Every operation is a subroutine taking an output `trend_result`.
Square brackets below denote optional arguments, not Fortran syntax.

| Procedure | Arguments after the name |
|---|---|
| `mk_test` | `x, result [, alternative, continuity]` |
| `itam` | `x, result [, conf_level]` |
| `sens_slope` | `x, result [, conf_level]` |
| `sea_sens_slope` | `x, period, result` |
| `smk_test` | `x, period, result [, alternative, continuity]` |
| `csmk_test` | `x, period, result [, alternative]` |
| `mult_mk_test` | `x, result [, alternative]` |
| `partial_mk_test` | `x, y, result [, alternative]` |
| `partial_cor_trend_test` | `x, z, result [, method]` |
| `bartels_test`, `cs_test`, `wm_test`, `ww_test`, `pettitt_test` | `x, result` |
| `rrod_test` | `x, y, result [, alternative]` |
| `lanzante_test` | `x, result [, method]` |
| `br_test`, `bu_test`, `snh_test` | `x, result [, m, seed, null_samples]` |

- `alternative`: exact strings `'two.sided'` (default), `'greater'`, `'less'`.
- `continuity`: default `.true.` for MK/seasonal MK. Multisite and correlated
  seasonal MK follow upstream: correction only with at most ten rows/cycles.
- `conf_level`: default `0.95_dp`, strictly between zero and one.
- `method` for partial correlation: `'pearson'` (default) or `'spearman'`.
- `method` for Lanzante: `'wilcox.test'` (default) or `'rrod.test'`.
- Multisite `x` has observations in rows and sites in columns, at least 3 by 2.
- Seasonal `x` is a consecutive vector with at least two seasons; `period`
  supplies seasons per cycle. SMK allows an incomplete terminal cycle with
  at least two observations per season. Correlated SMK requires complete cycles
  with at least three observations per season. Seasonal Sen requires complete
  cycles with at least two observations per season. First observation defines
  the first reported season; R start/end/calendar labels are not constructed.
- Ordinary Sen slopes are per observation interval; seasonal slopes are per
  complete cycle. Ordinary Sen accepts 3–10000 observations. Seasonal Sen limits
  each season to 10000 values and the combined slope array to 50 million pairs.
  These explicit bounds prevent integer-overflow-sized allocations.

### Result fields

`statistic` and `p_value` carry the principal test results. `score`, `variance`,
`tau`, `estimate`, `standard_error`, `correlation`, and `conf_int(2)` carry corresponding estimates
where applicable. Undefined scalar fields are quiet NaNs after each call.
`n` is retained observations (rows for multisite/correlated seasonal MK).

`itam` returns its slope in `estimate`, slope SD in `standard_error`, and
sorted first/last halves in `sorted_halves(:,1:2)`. With odd n the midpoint is
excluded from the halves but retained in the full-series SD and n denominator,
matching upstream. A zero SD with nonzero slope yields infinite z and p=0.

`smk_test` allocates `season_score`, `season_variance`, `season_z`, `season_p`,
and `season_tau`. Multisite/correlated seasonal tests allocate `covariance`;
partial correlation allocates the correlation matrix in `(time,x,z)` order in
that same field. Partial-correlation Student-t degrees of freedom are `n-2`,
matching this upstream implementation.

Change-point tests allocate `data` (the rank/cumulative/SNHT process) and
`change_points` (all maximizing one-based split indices). For br/bu, `data` is
always cumulative centered observations divided by **sample** standard deviation;
this matches upstream numeric-vector input. Upstream's inconsistent `ts` branch,
which returns unscaled cumulative sums for these two functions, is not reproduced.
Lanzante uses a unique interior split; tied maxima return `trend_ambiguous` rather
than reproducing R's vector-condition error. Its statistic is Wilcoxon W or RROD z.
`exact` identifies exact Wilcoxon inference. `simulations` reports null draws used.

`ww_test%variance` is the null variance after centering data on their first observation and
dividing by the maximum absolute centered value; this added diagnostic is scale-normalized to avoid overflow.
The reported z and p-value remain invariant to these transformations. This
corrects upstream raw-moment cancellation for large offsets: compare with R
`ww.test(x-x[1])` in that regime. For three observations the cyclic-product
statistic is invariant under permutation, so its zero variance is diagnosed.
Change-point statistics likewise center and scale before computing standardized
cumulative sums to improve stability for large offsets and magnitudes. ITAM
uncertainty and all Pearson correlations use normalized intermediates to avoid
overflow while retaining estimates and confidence bounds in original units.

### Status and special values

Check `status` after every call:

| Constant | Meaning |
|---|---|
| `trend_ok = 0` | Requested computation succeeded |
| `trend_invalid = 1` | Invalid options, dimensions, confidence level, or simulation count |
| `trend_insufficient = 2` | Too few observations |
| `trend_nonfinite = 3` | Nonfinite input in an API requiring finite data |
| `trend_degenerate = 4` | Zero/invalid variance or undefined correlation |
| `trend_ambiguous = 5` | Lanzante has multiple maximizing candidate splits |
| `trend_limit = 6` | Slope allocation safety limit exceeded |

MK, slope, multivariate, partial, Bartels, Pettitt, and simulation APIs reject
NaNs/infinities. `cs_test`, `wm_test`, and `ww_test` omit NaNs but reject retained
infinities. `rrod_test` independently omits all nonfinite observations in each
sample, exactly as upstream. A separated RROD sample legitimately produces
infinite z and p equal to zero or one; that is a successful result. Constant
Pettitt data return p=1 and every maximizing index. Constant MK/Bartels or
undefined correlations return a degeneracy status rather than relying on R
warnings or propagating accidental arithmetic exceptions. Seasonal constant
components have NaN season z/p/tau and zero score/variance; nonconstant seasons
can still yield a valid aggregate test.

An out-of-range Sen confidence rank returns NaN in that bound, rather than
reproducing R's zero/negative-index and variable-length vector behavior.
For a constant series its slope may be available with a degeneracy status.
Vector methods impose explicit minimum sizes where upstream R indexing would
otherwise misbehave (for example, WW requires at least three).

## Monte Carlo reproducibility

`br_test`, `bu_test`, and `snh_test` implement the actual upstream Monte Carlo
normal-null test, including **strict** exceedance `simulated > observed` and the
unadjusted count/m p-value. Default `m` is 20000. There is no plus-one correction.
As in upstream, the observed statistic uses sample standard deviation (divisor
n-1), but simulation statistics use population standard deviation (divisor n).

The local deterministic RNG is Park-Miller with Box-Muller normals. An optional
integer `seed` defaults to 104729; the same seed/data/options give the same result
(up to platform math-library rounding). It does not read or alter Fortran's global
RNG. The recurrence is not R's Mersenne-Twister/inversion generator; equal seed
numbers do **not** imply equal R and Fortran draws.

For cross-language reproducibility, pass `null_samples(n,m)` containing normal
samples, one realization per column. The column count becomes `m` unless an
explicit matching `m` is supplied. Input samples must be finite and nonconstant.
The deterministic tests use actual R `rnorm` samples and match all three upstream
p-values by exact exceedance counts. Additional independent-generator 50000-draw
comparisons check distributional agreement within Monte Carlo uncertainty.

## Translation coverage

Package status: **substantial**. **19 of 19 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility. S3 objects, formula evaluation, R warning/error
text, plotting, printing, recycling, calendar metadata, and R RNG streams are
not reproduced. Each mapping below is substantial because the numerical
operation is implemented with the documented interface and edge differences.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bartels.test` | `trend_module` | `bartels_test` | substantial |
| `br.test` | `trend_module` | `br_test` | substantial |
| `bu.test` | `trend_module` | `bu_test` | substantial |
| `cs.test` | `trend_module` | `cs_test` | substantial |
| `csmk.test` | `trend_module` | `csmk_test` | substantial |
| `itam` | `trend_module` | `itam` | substantial |
| `lanzante.test` | `trend_module` | `lanzante_test` | substantial |
| `mk.test` | `trend_module` | `mk_test` | substantial |
| `mult.mk.test` | `trend_module` | `mult_mk_test` | substantial |
| `partial.cor.trend.test` | `trend_module` | `partial_cor_trend_test` | substantial |
| `partial.mk.test` | `trend_module` | `partial_mk_test` | substantial |
| `pettitt.test` | `trend_module` | `pettitt_test` | substantial |
| `rrod.test` | `trend_module` | `rrod_test` | substantial |
| `sea.sens.slope` | `trend_module` | `sea_sens_slope` | substantial |
| `sens.slope` | `trend_module` | `sens_slope` | substantial |
| `smk.test` | `trend_module` | `smk_test` | substantial |
| `snh.test` | `trend_module` | `snh_test` | substantial |
| `wm.test` | `trend_module` | `wm_test` | substantial |
| `ww.test` | `trend_module` | `ww_test` | substantial |

The pinned NAMESPACE explicitly exports these 19 functions. The registered
`rrod.test.default`/`rrod.test.formula` methods are not independently exported;
the local numerical method is counted once under `rrod.test`. Plot/summary
methods, datasets, and internal helpers are outside the denominator. See
`API_COVERAGE.md`, per-function manifest notes, and `NOTICE.md`.

## Verification and license

`fpm test` runs stored actual-R reference values for all 19 operations, tied and
untied branches, tails, normal/beta transitions, seasonal covariance, all Monte
Carlo kernels, deterministic RNG, elementary independent oracles, invariances,
and invalid/special-value handling. See `validation/README.md` for regeneration
and checked/release commands. No R installation is required for these tests.

GPL-3.0-only for this combined translation; see `LICENSE` and `NOTICE.md`. Upstream source attribution and
license headers are retained. No dependency or build product is bundled.
