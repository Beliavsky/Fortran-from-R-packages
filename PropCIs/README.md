# PropCIs: modern Fortran confidence intervals

Translation of Ralph Scherer's **PropCIs 0.3-0** computational R code. Provides
single-binomial, paired-proportion, independent-difference, risk-ratio, odds-ratio,
and Bayesian beta-posterior intervals. GPL version 2 is elected from upstream's
unversioned `GPL` license declaration; see LICENSE and NOTICE.md.

## Build and run

Extract or place `PropCIs` directly under the root of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
with the existing sibling packages in place, then run:

```text
cd PropCIs
fpm build
fpm test
fpm run --example
```

Requires FPM and a modern gfortran. Uses free-form Fortran and the shared `dp` kind
from `rfortran-core`. There are no system BLAS/LAPACK links and no platform-specific
code, shell calls, or external runtime requirements for the Fortran tests/examples.
The source is portable to Windows gfortran; the recorded validation was on Linux,
not a separately tested Windows installation. Do not enable `-ffast-math`, `-Ofast`,
`-ffinite-math-only`, or equivalent flags: NaN and infinity are part of the API.

Sibling FPM dependencies, never copied into this package:

- `rfortran-core`: kind, normal quantile, beta/special functions (MIT)
- `extraDistr`: existing beta quantile and binomial distribution helpers (GPL-2.0-only)
- `RcppNumerical`: existing adaptive Gauss-Kronrod finite/infinite integration (GPL-2.0-or-later)
- `RcppNumerical` declares the transitive `lbfgs` (GPL-2.0-or-later) and `lbfgsb3`
  (GPL-2.0-only) siblings. No optimization or linear algebra sources are vendored.

The `stats` sibling was inspected but is not a dependency: at the tested repository
revision its build references an unavailable `r_matrix_statistics` module. Existing
focused dependencies provide the compatible APIs used here.

## API

```fortran
use propcis, only: dp, ci_result, exactci, diffci_bayes, ci_ok
implicit none
type(ci_result) :: ci
ci = exactci(5, 20, 0.95_dp)
if (ci%status /= ci_ok) error stop 'interval failure'
print *, ci%lower, ci%upper
ci = diffci_bayes(5, 20, 10, 30, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
```

All intervals return `type(ci_result)`:

- `lower`, `upper`: ordered interval bounds; infinite ratio endpoints are meaningful
- `estimate`: upstream point estimate for `add4ci`, `addz2ci`, the two paired Wald
  functions, and `wald2ci`; quiet NaN when upstream provides no estimate
- `conf_level`: requested confidence/posterior mass
- `status`: `ci_ok=0`, `ci_invalid=1`, `ci_numerical=2`

Invalid interval inputs return NaN bounds and `ci_invalid`; numerical failure can
return provisional finite bounds or NaNs with `ci_numerical`. Always check status.
The scalar exports `acceptbin`, `limit`, and `z2stat` return NaN for invalid inputs
or exhausted searches. The API does not print warnings or terminate on invalid data.

Fortran uses integer counts. Frequentist trial counts must be positive; paired
counts must be nonnegative with `b+c <= n`. Bayesian counts permit zero trials for
prior-only inference. Confidence levels are finite and strictly between zero and
one; all four beta prior shapes are positive and finite. The Bayesian argument
order is `(x1,n1,x2,n2,a,b,c,d,conf_level)`, with optional level defaulting to 0.95.
No simulation count or RNG argument is needed. `wald2ci(...,adjust)` accepts exactly
`"AC"` (default) or `"Wald"`. `blakerci(...,tolerance)` defaults to `1e-5` and accepts
`1e-9 <= tolerance <= 1`; a tiny step can be expensive.

`bayes_cdf(distribution,threshold,shapes,status,upper_tail,tolerance)` additionally
exposes the posterior CDF for independent verification. `shapes(4)` are *posterior*
shapes, not prior shapes; distribution is `bayes_difference`, `bayes_risk_ratio`, or
`bayes_odds_ratio`. The optional upper-tail switch computes survival probability
directly. Absolute quadrature tolerance defaults to `1e-10`.

### Numerical behavior and intentional differences

- `midpci` deliberately retains the R grid of 2,000 probabilities from 0.0001 to
  0.9996, rather than silently returning continuous roots. `blakerci` retains its
  requested additive grid; `orscoreci` retains the original multiplicative 1.001
  grid. `scoreci` retains four-decimal output, using `ANINT` rather than R's
  ties-to-even convention at exact decimal ties.
- `z2stat`/`diffscoreci` use a stable constrained-likelihood derivative root instead
  of the fragile cubic branch. `riskscoreci` brackets the same score equation;
  interior answers agree with the R cubic, while boundary answers can differ by
  the original 1.0001 grid spacing. Searches are bounded and failures explicit.
- Normal critical values are computed from the small tail, avoiding rounding a
  near-one confidence level to a probability of exactly one.
- Interior `exactci`, `blakerci`, and interior `midpci` mark confidence tails with
  `1-conf_level < 1e-10` as `ci_numerical`: beta-quantile resolution and binomial
  cancellation cannot establish accurate tail probabilities there. The closed-form
  zero/all-success `exactci` and `midpci` branches remain available.
- Bayesian equal-tail intervals integrate the exact independent beta-posterior
  target and invert it deterministically. The upstream beta/F Monte Carlo samples
  only initialize `optim`; those random starts, `nsim`, RNG state, and optimizer
  trajectories are intentionally omitted. Results are reproducible and need no seed.
  The odds-ratio integral is the equivalent beta-logit convolution, not a sampled
  F approximation. Standardized beta-logit coordinates regularize endpoint
  singularities and integrate over the full real line using the shared integrator.
- Bayesian inversion uses absolute quadrature tolerance at most `1e-11`, relative
  tolerance `1e-10`, and endpoint/log-ratio bisection width `2e-10`. Quadrature is
  limited to 1,000 subintervals. Shapes outside `[1e-6,1e6]`, requested excluded mass
  below `1e-8`, and ratio roots beyond log magnitude 700 return `ci_numerical`.
  These are numerical resource limits, not claims that inference is undefined.
  Strong endpoint singularities (for example all shapes 0.1) can also fail the
  quadrature tolerance and are reported, never silently accepted.
- `diffci_bayes_hpd` computes an approximate **shortest connected credible interval**
  by scanning 21 tail allocations and refining the best neighborhood with golden
  section search. A disconnected HPD set for multimodal/U-shaped posteriors is not
  returned. The scan is not a proof of the global minimum for arbitrary multimodal
  cases. For identical log-concave posteriors, the symmetric equal-tail interval is
  used exactly. On other cases, endpoint locations may vary by a few `1e-5` when
  the width objective is very flat, while posterior mass remains accurate.

### Verified upstream defects corrected

1. `orci.bayes` references `tailci` before assignment when `x2 == n2`. Direct posterior
   inversion handles that branch and retains asymmetric prior shapes correctly.
2. `diffci.bayes.hpd` clips the lower bound when its upper bound exceeds one and
   complements observations without swapping asymmetric beta prior shapes. Direct
   bounded posterior inference avoids both defects.
3. `wald2ci` passes its two-element default adjustment to `switch`, which errors
   without an explicit adjustment. Fortran selects `"AC"` by default.
4. `oddsratioci.mp` may produce huge finite or negative odds at a mathematically
   infinite endpoint due to roundoff. A zero first discordant count explicitly
   gives an infinite upper bound; a zero second count gives a zero lower bound.
5. The exported OR `limit` helper divides by zero at an initial ratio of one.
   Its equivalent stable quadratic formula is used instead.

R formula/S3/`htest` object construction, vector recycling, warnings, and optional
R interfaces are omitted. Scalar procedures can be called in user loops. Most
frequentist procedures are pure elemental and validation procedures are pure; beta quantile wrappers and Bayesian
procedures inherit impure interfaces from the reused dependency APIs.

## Translation coverage

**Package status: substantial. 21 of 21 (100%).** This fraction measures mapped
exported computational R functions, not complete R compatibility. The public
`acceptbin`, `limit`, and `z2stat` helpers count because the upstream namespace
exports their names. Every mapping is substantial with differences above and in
API_COVERAGE.md. `fpm.toml` records the same denominator and per-function mapping.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `acceptbin` | `propcis_frequentist` | `acceptbin` | substantial |
| `add4ci` | `propcis_frequentist` | `add4ci` | substantial |
| `addz2ci` | `propcis_frequentist` | `addz2ci` | substantial |
| `blakerci` | `propcis_frequentist` | `blakerci` | substantial |
| `diffci.bayes` | `propcis_bayesian` | `diffci_bayes` | substantial |
| `diffci.bayes.hpd` | `propcis_bayesian` | `diffci_bayes_hpd` | substantial |
| `diffpropci.Wald.mp` | `propcis_frequentist` | `diffpropci_wald_mp` | substantial |
| `diffpropci.mp` | `propcis_frequentist` | `diffpropci_mp` | substantial |
| `diffscoreci` | `propcis_frequentist` | `diffscoreci` | substantial |
| `exactci` | `propcis_frequentist` | `exactci` | substantial |
| `limit` | `propcis_frequentist` | `limit` | substantial |
| `midPci` | `propcis_frequentist` | `midpci` | substantial |
| `oddsratioci.mp` | `propcis_frequentist` | `oddsratioci_mp` | substantial |
| `orci.bayes` | `propcis_bayesian` | `orci_bayes` | substantial |
| `orscoreci` | `propcis_frequentist` | `orscoreci` | substantial |
| `riskscoreci` | `propcis_frequentist` | `riskscoreci` | substantial |
| `rrci.bayes` | `propcis_bayesian` | `rrci_bayes` | substantial |
| `scoreci` | `propcis_frequentist` | `scoreci` | substantial |
| `scoreci.mp` | `propcis_frequentist` | `scoreci_mp` | substantial |
| `wald2ci` | `propcis_frequentist` | `wald2ci` | substantial |
| `z2stat` | `propcis_frequentist` | `z2stat` | substantial |

## Verification

`fpm test` includes 255 assertions against actual PropCIs 0.3-0 outputs, 18 posterior
interval cases from independent R `integrate`/`uniroot` distribution definitions,
and deterministic domain, boundary, confidence-tail, known-defect, sample-exchange,
credible-mass, shortest-width, and analytic uniform-prior tests. Cases include zero
and all successes, highly unequal proportions, Jeffreys and skew fractional priors,
asymmetric priors, small tails, and finite/infinite ratio support boundaries.
The all-success R Bayesian-OR error is recorded and its corrected result tested
against the independent posterior definition rather than a fictitious R result.

Checked-in `test/*reference.tsv` contain the numeric oracles. Optional regeneration
requires R with the unmodified PropCIs 0.3-0 package and Python 3:

```text
Rscript tools/generate_reference.R
python tools/embed_reference.py
```

R is not required to build or run the deterministic Fortran tests. Reference
regeneration deliberately records upstream errors/warnings instead of hiding them.
See VALIDATION.md for the build/test/example configurations and scope.
