# Porting notes

## Architecture

The R package contains many nearly identical files. The Fortran translation
factors them into:

- `bgfd_core`: eight baseline distributions and the two Bell transforms
- `bgfd_distributions`: named wrappers corresponding to the R d/p/q/r/s/h API
- `bgfd_fit`: common MLE/goodness-of-fit engine plus the sixteen `m_*` wrappers
- `bgfd`: public facade module

All distribution parameters are required to be positive and all distributions
have support on `x >= 0`.

## Corrected Weibull quantiles

The upstream functions `qBellW`, `qCBellW`, `qBellEW`, and `qCBellEW` apply the
power `beta` after solving for `x^beta`. The inverse of the CDF requires the
power `1/beta`. The upstream expressions therefore fail the identity
`F(Q(p)) = p` whenever `beta /= 1`.

The Fortran routines use `1/beta`. Permanent tests check both the independent
reference values and CDF/quantile inversion for all sixteen families.

## Survival and hazard flag semantics

The upstream `s*` routines calculate a CDF or log-CDF according to `log.p` and
`lower.tail`, then always apply `1-cdf`. Consequently `log.p=TRUE` does not
return a log-survival probability. Similarly, the `h*` routines divide the
quantity returned by the density branch by `1-cdf`; when `log=TRUE` that means
a log-density is divided by a probability rather than returning a log-hazard.

The Fortran API uses conventional semantics:

- `s_*` returns survival probability; `log_value=.true.` returns log-survival.
- `h_*` returns the hazard; `log_value=.true.` returns log-hazard.
- `p_*` retains lower/upper-tail and log-probability controls.
- `q_*` retains lower/upper-tail and log-probability controls.

Default (non-log, lower-tail) density/CDF/survival/hazard behavior matches the
mathematical formulas used by the R package.

## Numerical stability

The Bell transforms are evaluated with small-argument `expm1`/`log1p` analogues
and log-domain complementary-Bell normalizing constants. This avoids avoidable
cancellation and postpones overflow for moderately large `lambda`.

Support boundaries are handled explicitly rather than relying on expressions
such as `0**negative` to signal a limiting density.

## Fitting

The upstream `m*` functions delegate to `AdequacyModel::goodness.fit`. The
Fortran port vendors the compatible AdequacyModel translation and preserves the
same optimizer choices and goodness-of-fit outputs.

The positive BGFD parameters are optimized as their logarithms. This preserves
the target likelihood while preventing optimizers from evaluating invalid
negative shape/scale parameters. Natural-scale standard errors are obtained by
the delta method. If the default BFGS optimizer stalls before satisfying its
finite-difference gradient criterion, the solution is refined from the same
point with Nelder-Mead.

The fit interface is module-state based while callbacks are active, matching the
callback limitations of the vendored numerical layer; concurrent fits from
multiple threads should therefore use external synchronization.
