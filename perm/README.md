# perm: modern Fortran permutation tests

A computational translation of Michael P. Fay's R package **perm 1.0-0.4**:
exact two-sample tests, complete enumeration, grouped-score network enumeration,
Monte Carlo two-sample/k-sample/trend tests, asymptotic approximations, and
p-value uncertainty intervals. **permKS is a k-sample test, not a
Kolmogorov–Smirnov test.**

## Build and run

Place `perm/` beside the repository's `rfortran-core/`, then run:

```sh
fpm build
fpm test
fpm run --example --all
```

Requires a modern gfortran and FPM. All maintained Fortran is free form, uses the
single `dp` kind imported from the sibling core, and has explicit interfaces. No
system BLAS/LAPACK, R installation, C compiler or separately installed numerical
library is needed for the Fortran build. The path dependency works on Windows
without platform-specific library links. R is needed only to regenerate reference
fixtures. Do not use `-ffast-math`, `-Ofast`, `-ffinite-math-only` or equivalent
flags: validation relies on IEEE finite/NaN semantics.

## Quick example

```fortran
use perm_api, only : dp, perm_result, perm_ts, perm_ok
type(perm_result) :: result
call perm_ts([1.0_dp, 1.0_dp, 2.0_dp], [3.0_dp, 5.0_dp, 8.0_dp], &
   result, method='exact.network')
if (result%status /= perm_ok) error stop 'Permutation test failed'
print *, result%p_value              ! 0.1
print *, result%estimate             ! -4: first mean minus second mean
```

All numerical procedures are `pure`; scalar numerical helpers are `pure elemental`
where possible. Caller-owned RNG state makes Monte Carlo calls pure as well.
Every dummy argument has its own declaration and trailing FORD documentation.

## API conventions

Import the public module `perm_api`. Tests return `type(perm_result)`:

- `status`: always check it before interpreting fields.
- `p_value`: selected probability; `p_twosided`, `p_absolute`, `p_less`,
  `p_greater`, `p_equal` expose the individual tails when meaningful.
- `statistic`: normal Z or chi-square for asymptotic tests; observed score sum
  for exact enumeration. MC statistics use internal power-of-two normalized
  scores/covariates; upstream R does not expose a numeric MC statistic.
- `estimate`: mean(x)-mean(y) for `perm_ts`, correlation for `perm_trend`.
- `df`: k-sample asymptotic degrees of freedom, otherwise zero.
- `confidence_interval(2)` and `nmc`: meaningful for Monte Carlo results only.
- `method`: `exact.network`, `exact.ce`, `exact.mc`, or `pclt`.

A probability/interval field not supplied by the selected method remains `-1`.
An invalid result is not a valid p-value; its fields are unspecified. A constant
sample/covariate makes asymptotic inference undefined and returns
`perm_degenerate`. Exact and MC tests can validly return one for all-tied data;
the undefined trend correlation is NaN.

`perm_ts(x,y,...)` takes the first sample as the positive/upper-tail sample.
Low-level `twosample_pclt`, `twosample_exact_ce` and
`twosample_exact_network` take exactly two integer labels and select the larger
label. `twosample_exact_mc` requires explicit zero/one labels, as its R kernel
does. K-sample labels can be any integers and are compressed into sorted groups.
Trend covariates are real arrays. Missing/NaN/infinite observations are rejected;
there is no implicit omission or recycling.

High-level tests have optional `alternative`, `exact`, `method`, `control` and/or
`rng` arguments as appropriate. `exact=-1`, `0`, `1` mean R's `NULL`, `FALSE`,
`TRUE`. An explicit method overrides automatic selection. A custom R method-rule
callback is replaced by calling your rule externally and passing the resulting
method name. Formula processing, R object construction and S3 presentation are
omitted.

## Exact algorithms, ties and tails

`choose_matrix` enumerates selection rows in R's order, including complement
ordering when m>n/2. `twosample_exact_ce` evaluates all such rows. A supplied
`control%cm` must be binary, have the correct column count and select the observed
number of observations in each row. As in R, supplied rows are taken literally;
completeness and uniqueness are the caller's responsibility.

`twosample_exact_network` compresses equal score values and enumerates feasible
occupancy paths. The probability of each path is proportional to the product of
its binomial multiplicities. Log weights avoid factorial/binomial overflow. No
score lattice or integer rounding approximation is used. Consequently arbitrary
finite real scores and ties are supported within the stated resource bounds.

`calc_pvals` normalizes nonnegative weights, centers all statistics at their
weighted mean, and compares centered statistics rounded to `digits` significant
decimal digits (default 12). Decimal half ties use ties-to-even. A decimal carry
is canonicalized to avoid unequal floating representations of equal rounded
magnitudes. Extremely small/large statistics are rescaled before rounding; subnormal-scale
statistics are also normalized before weighted centering to avoid underflow of
weights times scores. Like
R, outcomes exactly at a rounding boundary can depend on floating-point arithmetic
and summation order; use appropriate score precision and `digits` for the data.
The unused R `denom` argument is omitted: R itself overwrites it with sum(weights).

The two-sided conventions are different for asymmetric permutation distributions:

- `central`: min(1, 2*min(P(T<=observed), P(T>=observed))).
- `abs`: P(abs(T-mean(T))>=abs(observed-mean(T))).

`alternative='two.sided'` uses `control%tsmethod`, which defaults to `central`.
`two.sidedAbs` is also accepted. `less`/`greater` include equality. An already
two-sided statistic, such as k-sample between-group variation, uses its upper tail.

## Monte Carlo and uncertainty

`perm_control` builds `perm_control_type`. High-level calls default to 999 random
permutations; the standalone two-sample and k-sample MC kernels default to 9,999,
and the standalone trend kernel defaults to 999, matching upstream defaults.
The observed statistic is always added: p=(tail_hits+1)/(nmc+1). Centering for the
absolute convention uses the simulated mean, as in R.

Sampling is Fisher–Yates with unbiased bounded-integer rejection draws from a
local Park–Miller generator (multiplier 48271, modulus 2147483647). Seeds are
normalized by 1+modulo(seed,2147483646). The default is 1234321. It deliberately
**does not reproduce R's RNG sequences** and never changes Fortran's global RNG.
The same inputs and seed reproduce results across supported builds. With
`set_seed=.true.` each call resets its local stream. With `set_seed=.false.`, pass
an explicit `type(perm_rng)` whose state is in 1..2147483646; it is updated after a
successful call. There is no implicit process-global continuation.

Each simulation runs its full fixed `nmc`. There is **no sequential stopping** or
claim of a confidence sequence. `p_ci` computes a fixed-count equal-tailed
Clopper–Pearson interval from the raw tail hits, before the add-one correction.
Central two-sided intervals double the selected smaller-tail interval and cap at
one. As upstream, absolute-tail indicators use the simulated center; the interval
is a Monte Carlo diagnostic, not an anytime-valid guarantee for data-dependent
centering or repeated optional stopping.

Two numerical repairs are deliberate:

1. Both doubled central interval endpoints are capped at one. R can leave the
   lower endpoint above one for all-tied draws, yielding a reversed interval.
2. Small beta tails use stable symmetry, analytic shape-one formulas, and shared
   F-quantile/beta-CDF refinement. At confidence immediately below one, this
   avoids rounding an upper-tail probability to exactly one as R can do.

MC inputs are power-of-two scaled before products/quadratic statistics. This
preserves the permutation statistic's ordering while preventing tiny scores from
underflowing to an all-tied distribution. R seed handling, warnings and global
state are not emulated.

## Status codes and practical bounds

| Status | Meaning |
|---|---|
| `perm_ok` = 0 | Computation completed |
| `perm_invalid` = 1 | Invalid shape, control, category, or nonfinite input |
| `perm_limit` = 2 | A declared resource limit was exceeded |
| `perm_allocation` = 3 | An explicitly checked allocation failed |
| `perm_degenerate` = 4 | Zero variance makes asymptotic inference undefined |
| `perm_numeric` = 5 | Magnitudes exceed the safe numerical range |

The default limits are 1,000,000 exact paths, 10,000,000 principal table elements,
and 100,000,000 estimated enumeration/MC operations. `max_permutations`,
`max_cells` and `max_work` in controls can change these caps. Sorting/grouping
scratch arrays are O(n) in addition to principal tables. The hard limits are
100,000 observations, 10,000,000 MC trials, 512 distinct-score groups in the
recursive network traversal, and 1..15 significant digits. The singleton-group
network shortcut also honors path/cell/work caps. Complete enumeration can handle
more than 512 distinct scores if its separate combinatorial limits allow it.
Counts are bounded before conversion/allocation, and bounded-combination arithmetic
uses int64. Controls do not silently switch an explicitly requested exact method
to MC or asymptotic inference.

Safe finite-input magnitude for score arrays is sqrt(huge(dp))/(8*n). This
conservative documented boundary prevents arithmetic overflow in all kernels.
The weighted-tail utility accepts a larger range, up to huge(dp)/8 for statistics.
Confidence may be any finite value in [0,1]; at one the interval is [0,1].

Automatic TS selection retains upstream distinct-score thresholds, except that a
group with one distinct score uses the first threshold instead of upstream's
invalid zero index. Resource-limit failures remain explicit even if an automatic
rule selected an exact method.

## Validation

`fpm test` exercises:

- Numerical fixtures generated by the official R package: weighted tails, both
  exact algorithms, all asymptotic kernels, both two-sided MC conventions,
  pCI boundaries and large counts, selection-row ordering and edge indices.
- Independent enumeration of every nontrivial group assignment for n=3..8,
  across three tied-score families; complete-enumeration/network agreement.
- 100 deterministic random real-score inputs, scale/location and sample-swap
  identities; all tails checked separately.
- 100 fixed-seed 1,999-draw simulations against exact probabilities and 99%
  interval coverage; independent six-factorial k-sample/trend calibration.
- Reproducible/continued RNG state, tiny-scale MC invariance, all-ties intervals,
  invalid inputs, rounding, degeneracy, combinatorial and operation limits.

Regenerate the stored fixture with an installed official R package:

```sh
Rscript reference/generate_reference.R
```

The example covers exact, k-sample MC, trend asymptotic and explicit-stream tests.
Default, checked (`-fcheck=all`) and release builds/tests/examples are required for
release validation. See `NOTICE.md` for exact source provenance and numerical
reuse decisions. Full upstream sources and licenses are retained.

## Translation coverage

Package status: **substantial**. **20 of 20 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility. Numerical validation, resource bounds, RNG behavior,
formula/S3 omissions and the deliberate repairs above remain material differences.
`permControl` and `getcnt` are included for computational validation and network
index generation. Registered but unexported S3 default/formula methods are not
separately counted; unexported `print.mchtest` is presentation-only.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `chooseMatrix` | `perm_api` | `choose_matrix` | complete |
| `permTS` | `perm_api` | `perm_ts` | substantial |
| `permKS` | `perm_api` | `perm_ks` | substantial |
| `permTREND` | `perm_api` | `perm_trend` | substantial |
| `getcnt` | `perm_api` | `getcnt` | complete |
| `methodRuleTS1` | `perm_api` | `method_rule_ts1` | substantial |
| `methodRuleKS1` | `perm_api` | `method_rule_ks1` | substantial |
| `methodRuleTREND1` | `perm_api` | `method_rule_trend1` | substantial |
| `twosample.exact.ce` | `perm_api` | `twosample_exact_ce` | substantial |
| `twosample.exact.mc` | `perm_api` | `twosample_exact_mc` | substantial |
| `twosample.exact.network` | `perm_api` | `twosample_exact_network` | substantial |
| `twosample.pclt` | `perm_api` | `twosample_pclt` | substantial |
| `ksample.pclt` | `perm_api` | `ksample_pclt` | substantial |
| `ksample.exact.mc` | `perm_api` | `ksample_exact_mc` | substantial |
| `trend.pclt` | `perm_api` | `trend_pclt` | substantial |
| `trend.exact.mc` | `perm_api` | `trend_exact_mc` | substantial |
| `calcPvals` | `perm_api` | `calc_pvals` | substantial |
| `permControl` | `perm_api` | `perm_control` | substantial |
| `calcPvalsMC` | `perm_api` | `calc_pvals_mc` | substantial |
| `pCI` | `perm_api` | `p_ci` | substantial |

See `API_COVERAGE.md` and the exact per-export mappings in `fpm.toml`.

## License

GPL-2.0-only OR GPL-3.0-only, preserving the upstream choice. See `LICENSE`,
`LICENSE-GPL-2`, `LICENSE-GPL-3` and `NOTICE.md`. The sibling core is separately
MIT-licensed and is not bundled or modified here.
