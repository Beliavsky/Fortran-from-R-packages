# Porting notes

## Scope

This is a translation of the computational code in PDQutils 0.1.6.  All nine
exports in the supplied `NAMESPACE` are represented.  Plotting used only in
examples/vignettes was not translated.

## Removed runtime dependencies

The R package imports `orthopolynom`, `moments`, and several `stats` routines.
The Fortran translation is standalone:

- probabilists' Hermite polynomials use their three-term recurrence;
- generalized Laguerre polynomials use their standard three-term recurrence;
- Jacobi polynomials use their standard three-term recurrence;
- normal PDF/CDF/quantile routines are provided internally;
- regularized incomplete gamma and beta routines provide gamma/beta CDFs;
- moment/cumulant conversion is translated directly.

This avoids imposing the licensing/dependency chain of a separate
`orthopolynom-fortran` build on PDQutils.

## Generalized Gram-Charlier beta CDF correction

In `R/gram_charlier.r`, the `intpoly` closure for the beta basis constructs a
Jacobi polynomial with parameters `(alpha+1,beta+1)` but evaluates the beta
weight with the corresponding shape parameters reversed.  The parent-beta
regression test cannot expose this because every correction coefficient is
zero when the target distribution equals the parent.

The Fortran port uses the mathematically consistent identity

```text
integral w(x) P_n^(alpha,beta)(x) dx
 = -B(alpha+2,beta+2)/(n B(alpha+1,beta+1))
   * beta_pdf((x+1)/2; beta+2, alpha+2)
   * P_(n-1)^(alpha+1,beta+1)(x).
```

A non-parent beta expansion in `test_gca.f90` checks the resulting CDF against
independent direct integration; the agreement is at floating-point roundoff.

## Arcsine endpoint normalization

The generic closed formula for the beta/Jacobi squared norm has a removable
singularity at degree zero for the arcsine parent (`alpha=beta=-1/2`).  The
zeroth norm is exactly one because the parent weight is a normalized density.
The Fortran implementation handles degree zero explicitly.  This makes the
arcsine basis usable rather than propagating an indeterminate `Inf-Inf` log
expression.

## Cornish-Fisher endpoints

Upstream `qapx_cf` passes probability endpoints through `qnorm` and then into
high-order polynomial corrections; nonzero higher cumulants can turn an
infinite normal quantile into `NaN`.  The Fortran routine returns the requested
support endpoint directly for probabilities 0 or 1 (and their log-probability
equivalents).  Interior probabilities follow AS269.

## Approximation semantics

As upstream documents, the Edgeworth and generalized Gram-Charlier CDFs need
not be monotone, and Cornish-Fisher quantiles need not be monotone.  For
compatibility, densities are clipped below at zero and CDFs are clipped to the
unit interval.

## Validation

The test suite contains:

- moment/cumulant round trips;
- independent AS269/Cornish-Fisher reference values;
- independent Edgeworth values for chi-square cumulants;
- exact-parent normal, gamma, beta, arcsine, and Wigner GCA cases;
- a non-parent beta GCA case;
- support/tail checks;
- vector interface checks;
- a Monte Carlo check of `rapx_cf` for a normal target.
