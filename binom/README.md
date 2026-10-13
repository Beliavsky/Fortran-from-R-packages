# binom: binomial inference in modern Fortran

Translation of CRAN **binom 1.1-2** by Sundar Dorai-Raj and Arianna D. Cascone.
Eleven confidence-interval methods, exact coverage/expected-length enumeration,
power and sample-size design, reproducible Monte Carlo summaries, and
integrated-coverage L-BFGS-B optimization. See [NOTICE.md](NOTICE.md) for source
checksums, dependency inspection and licenses.

## Build and run

Place this directory at the root of Beliavsky/Fortran-from-R-packages alongside
`rfortran-core`, `rfortran-linalg` and `lbfgsb3`. These are sibling FPM path
dependencies; none is included in this package. Use:

```
fpm build
fpm test
fpm run --example confidence_intervals
```

The source is intended to be portable across Windows and Linux with gfortran.
Validation for this translation was run on Linux; no native Windows run was
performed. No separately
installed BLAS, LAPACK, R, Python or C library is needed. The existing
rfortran-linalg dependency retrieves its tested portable fortran-lapack revision.
There are no system-library links. Its first source build is large and can take
several minutes. R is only needed to regenerate optional reference fixtures.

All maintained real code uses the single `dp` kind imported from `r_kinds`.
Public numerical routines are pure/elemental wherever their interfaces permit;
the shared optimization solver is the necessary impure exception. Do not compile
with `-ffast-math`, `-Ofast`, `-ffinite-math-only` or equivalent assumptions:
NaN/Inf validation is part of the documented behavior.

## Quick use

```
use binom
implicit none
type(binom_interval) :: ci
type(binom_design_result) :: design

ci = binom_confint(3.0_dp, 10.0_dp, 'wilson', 0.95_dp)
if (ci%status /= binom_ok) error stop 'Invalid interval'
print *, ci%mean, ci%lower, ci%upper

design = cloglog_sample_size(p_alt=0.7_dp, p=0.4_dp, recompute_power=.true.)
print *, design%n, design%power
```

`use binom` exposes the complete public API. Direct module use is also supported.
`binom_confint(x,n,method,conf_level)` returns one `binom_interval`; the default
method is Wilson. `binom_confint_all(x,n,conf_level)` returns 11 intervals in the
order of the public `binom_methods` array. For arrays of observations use the
elemental selected-method procedures; there is no implicit R recycling.
Counts and trial counts in interval routines are real(dp), permitting the
fractional counts used by upstream power calculations. Enumeration and
simulation routines take integer trial counts.

The eight direct wrappers `binom_wilson`, `binom_agresti_coull`, `binom_exact`,
`binom_asymp`, `binom_prop_test`, `binom_cloglog`, `binom_logit` and `binom_probit`
have the same `(x,n,conf_level)` arguments. Select the direct `binom_bayes`,
`binom_lrt` or `binom_profile` routine for method-specific options.

### Results and validation

- `binom_interval`: `mean`, `lower`, `upper`, `shape1`, `shape2`, `sig`, `status`.
  The shape and excluded posterior mass fields apply to Bayes only. Bayes
  `mean` is the posterior mean; other methods return the observed proportion.
- Interval status: `binom_ok=0`, `binom_domain=1`, `binom_convergence=2`,
  `binom_bad_method=3`, `binom_resource=4`. Invalid interval inputs give NaN
  endpoints. Incomplete iterative fits retain their last computed values with nonzero
  status. Always check status before using a structured result.
- All confidence probabilities must be strictly inside (0,1), n must be positive,
  x must be finite and in [0,n]. Invalid inputs are returned, never silently
  removed. Named methods are case-sensitive exact strings; `ac` and `asymp`
  aliases are supported. No R partial matching, data frames, plotting, S3,
  warnings, callbacks or arbitrary ellipsis arguments are reproduced.
- `binom_coverage` and `binom_length` return NaN on bad inputs or failed intervals.
  Their inclusive endpoints and unclipped interval widths follow upstream.
  Agresti-Coull and asymptotic limits intentionally may lie outside [0,1].
  Wilson limits at x=0 and x=n are set to their exact mathematical 0/1 values.
  This intentionally fixes upstream ulp-sized endpoint errors that can make
  R report zero coverage at true p=0 or p=1 for some trial counts.
- Enumeration and uniform simulation are capped at 10,000 trials. Profile
  accepts at most 4,096 requested steps per direction and 512 retained knots;
  exceeding a cap returns a resource status. These caps are explicit defenses
  against unbounded CPU/memory requests, not statistical limitations.

### Bayesian, likelihood and profile intervals

`binom_bayes` supports `interval_type='highest'` (default) or `'central'`, positive
`prior_shape1`/`prior_shape2` (both default 0.5), `tol` and `maxit`. As upstream,
observed endpoints get one-sided posterior intervals. Interior HPD intervals
require posterior shapes >1, always true for positive priors and integer
interior binomial observations. HPD uses equal log density and fixed posterior
mass rather than nested C density-height roots. The default relative
tail-allocation tolerance is 1e-12 and the iteration cap is 200. This is more
accurate than R's default absolute C tolerance, so small differences are expected.

`binom_lrt` supports optional `bayes` (default false), `conf_adj` (default false),
custom positive prior shapes, `tol` and `maxit` (default 1000, maximum 1000000).
Unrepresentable adjusted pseudocount sums return resource status before addition.
It uses the upstream Brent iteration and preserves a subtle R default: the
boundary/MLE-neighborhood tolerance is sqrt(epsilon), while the root solver's
omitted tolerance is epsilon**0.25, because R's local default is not forwarded
to uniroot. An explicit tol sets both values, matching R's ellipsis override.
Use a tight explicit tolerance when accurate likelihood-root residuals matter.
Boundary clipping and conf_adj tail changes follow the R branch logic. Stable
log-likelihood differences avoid subtracting the common binomial coefficient.

`binom_profile` is the actual sampled likelihood-profile method, not an alias
for LRT. It retains logit stepping, Jeffreys edge adjustment (default true),
FMM cubic interpolation on a 3-times-knot grid, and inverse linear interpolation.
`maxsteps` defaults to 50; omitted `del` is recomputed as zmax/5. A missing target
bracket returns convergence status and a NaN bound. Natural-cubic interpolation
from the existing splines package would change the algorithm and is not used.
In R, the lazy default del is evaluated once using the first vector observation.
Coverage, length and simulation preserve that shared step for their collections;
separate scalar/elemental profile calls recompute it per observation. Optimizing
profile intervals uses each count's scalar confidence level and does not reproduce
R's accidental vector recycling inside its profile loop.

### Coverage, power and design

`binom_coverage(p,n,method,conf_level)` sums the binomial masses whose intervals
contain p. `binom_length` sums mass times interval width. Both are elemental;
select multiple methods or parameter combinations with explicit arrays/loops.

`binom_power(p_alt,n,p,alpha,phi,alternative,method)` supports cloglog, logit,
probit, asymp, lrt and exact, with two.sided/greater/less alternatives. Defaults
match R (n=100,p=.5,alpha=.05,phi=1,cloglog,two.sided). Transformed-normal methods
with nonlinear links require p and p_alt strictly in (0,1). Asymptotic power
supports endpoint probabilities through degenerate-normal limits. The exact/lrt methods retain upstream's
beta-tail formula: despite the name “exact”, this is not a conventional integer
rejection-region power calculation. Degenerate beta shapes are supported at
p_alt endpoints. Upstream power can slightly exceed one there; it is not clipped.
Scalar/elemental evaluation avoids the R function's vector indexing/recycling bug.

`cloglog_sample_size` returns a `binom_design_result` with `p_null`, `p_alt`,
`delta`, `alpha`, `power`, `n`, `phi`, `status`:

- Supply p_alt and omit n to compute sample size (power defaults to .8).
- Supply n and p_alt to compute attained power.
- Supply n and omit p_alt to find a detectable alternative for the target power.
  Two-sided inversion searches above p; less searches below p. As in R, a
  requested target can fail to be bracketed by the endpoint probabilities.
  This inverse uses 100 bisection steps, rather than R uniroot’s default
  epsilon**0.25 stopping tolerance. The resulting probability can differ by
  a few 1e-5 while satisfying the requested power equation more accurately.
- `exact_n=.true.` returns a fractional size, at least one; otherwise the size
  rounds upward using real arithmetic without default-integer overflow.
- `recompute_power=.true.` updates power after rounding a computed size.

### Reproducible simulation

`binom_sim(uniforms,n,p,method,conf_level)` transforms caller-supplied uniforms
in [0,1) into binomial counts using beta-CDF inversion. `binom_sim_counts`
accepts counts directly. Both return a `binom_sim_result` containing the capture
count, replicate count, capture-proportion interval under the same method,
mean interval width, and status. Supplying random inputs separates numerical
results from compiler-dependent global RNG state. R's rbinom stream, default
M=200 and implicit draws are deliberately not emulated.

### Coverage optimization

`call binom_optim(n,fit,conf_level,method,k,p0,transform,start,control)` uses the
existing sibling lbfgsb3 solver. Defaults follow R's loss and symmetric expansion:
method lrt, k=n/2+1, p0=0, transform=true, starting confidence .95. The unused
bounds computed by the R source are not imposed. Untransformed optimization
retains R's closed-form interval algebra outside ordinary (0,1) confidence,
including signed-normal/reversed intervals and temporary NaN endpoint/count
removal. Public interval entry points retain their strict domain. A nonstandard
optimized level is an algebraic result, not a conventional confidence interval.
Trials with no finite interval are rejected instead of reproducing R's spurious
zero loss. A segment with no finite selected counts fails as in R; a solver stop
with invalid final bounds also returns status 53.
`start` is recycled across k.
`binom_optim_control` provides memory=5, max_iterations=100, max_evaluations=10000,
factr=1e7, pgtol=0 and finite_difference_step=1e-3. History memory must be in
1:1000 (`binom_optim_max_memory`); the finite-difference step must be positive
and at most huge(dp)/2. Violations return status 52. Iteration and total objective
caps are strict. Supported n is 1:200.

`binom_optim_result` supplies `par(k)`, `conf_level(0:n)`, `confint(0:n)`,
`initial_value`, `value`, convergence code, counters and message. Convergence
0 is success; 1 is a resource cap; 51 a solver warning; 52 invalid input;
53 a nonfinite loss/gradient. Check convergence even if a usable last iterate
was retained. L-BFGS-B versions and floating-point differences can change stopping codes at
nonsmooth interval crossings. The n=20 raw Agresti-Coull reference reaches the
same optimum but may return warning 51 instead of R's success code. This warning
is preserved rather than silently reclassified as convergence.
Pure helpers `binom_optim_objective` and `integrate_binomial_error` expose the
piecewise exact-beta-integral loss. R's optional polynom representation computes
the same integral and is not a dependency. General optim parscale/fnscale,
analytic user gradients, arbitrary R callback methods and plotting are omitted.

## Validation

`tools/generate_references.R` regenerates checked-in deterministic CSV fixtures
using the actual installed binom 1.1-2 (R 4.5.0 for this port). The suite includes
660 interval rows across all methods and multiple confidence levels, 66 exact
coverage/length comparisons, 72 powers, central/custom-prior HPD intervals,
adjusted LRT, custom profile steps, all sample-size branches, and deterministic
simulation comparisons. Additional tests check beta inverse identities/tails,
posterior mass/equal endpoint density, likelihood-root residuals, binomial mass
normalization, exact coverage, invalid inputs, NaN/Inf and resource caps.
Optimizer tests include actual-R baseline fits/objectives, independent integral
identities, symmetry, start handling and strict resource limits.

R's C HPD stopping tolerance justifies endpoint/power fixture tolerances
up to 3e-7; the Brent defaults are matched, and profile-specific fixtures use
3e-10. Independent likelihood-root checks request an explicit tight tol. Other identity checks are
substantially tighter. `fpm test` needs neither R nor a network after dependencies
are available. The reference-generation script is optional.

## Translation coverage

Status: **substantial**. **18 of 18 (100%)** exported computational R functions
are mapped. This is a function-mapping fraction, not a claim of complete R
compatibility. Three plotting/interactive exports are excluded from the
computational denominator. [API-COVERAGE.md](API-COVERAGE.md) records all mapping
limitations; this table agrees with `fpm.toml`.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `binom.asymp` | `binom_intervals` | `binom_asymp` | complete |
| `binom.agresti.coull` | `binom_intervals` | `binom_agresti_coull` | complete |
| `binom.exact` | `binom_intervals` | `binom_exact` | complete |
| `binom.wilson` | `binom_intervals` | `binom_wilson` | complete |
| `binom.prop.test` | `binom_intervals` | `binom_prop_test` | substantial |
| `binom.bayes` | `binom_intervals` | `binom_bayes` | substantial |
| `binom.cloglog` | `binom_intervals` | `binom_cloglog` | complete |
| `binom.probit` | `binom_intervals` | `binom_probit` | complete |
| `binom.profile` | `binom_profile_module` | `binom_profile` | substantial |
| `binom.confint` | `binom_intervals` | `binom_confint`, `binom_confint_all` | substantial |
| `binom.coverage` | `binom_design` | `binom_coverage` | substantial |
| `binom.lrt` | `binom_intervals` | `binom_lrt` | substantial |
| `binom.logit` | `binom_intervals` | `binom_logit` | complete |
| `cloglog.sample.size` | `binom_design` | `cloglog_sample_size` | substantial |
| `binom.power` | `binom_design` | `binom_power` | substantial |
| `binom.sim` | `binom_design` | `binom_sim`, `binom_sim_counts` | substantial |
| `binom.length` | `binom_design` | `binom_length` | substantial |
| `binom.optim` | `binom_optimization` | `binom_optim` | substantial |
