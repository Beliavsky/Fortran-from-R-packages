# Porting notes

## Scope

This project translates numerical/statistical code from evd 2.3-7.1 and does
not attempt to reproduce R's S3 classes, formula evaluation, plotting,
interactive diagnostics, profile plotting, print methods, or data-frame
presentation machinery.

The package's core distribution, dependence, simulation, stochastic-process,
nonparametric and stationary likelihood calculations are represented in
Fortran. R wrappers whose sole role is rearranging arguments or constructing
S3 return objects are replaced by ordinary typed Fortran APIs.

## `r_mod.f90`

The user-supplied MIT-licensed `r_mod.f90` is copied into `src/` and reused for
R-compatible numerical helpers, random generation and BFGS optimization.
Only evd-specific helpers were added.

## Bivariate POT likelihoods

The upstream C source writes separate closed-form censored likelihood routines
for each dependence family. The Fortran translation factors their common
mathematics through the homogeneous exponent representation

`V(x,y) = (x+y) A(x/(x+y))`.

For censored likelihoods, derivatives of `A` are evaluated numerically, while
the same marginal transformations and censoring contributions are retained.
For the logistic reference case, the translated value differs from the
upstream C formula by about 2.7e-9 on the included reference data.

For Poisson-process likelihoods the translated code uses the analytic angular
(`h`) functions already translated from evd. The logistic reference value
matches the upstream formula to displayed double precision.

As in upstream evd, asymmetric logistic, asymmetric negative-logistic and
asymmetric mixed models are not accepted by the Poisson likelihood.

## Multivariate asymmetric logistic density

Rather than copying a dimension-specific expanded expression, `dmvalog`
computes the exact mixed derivative using set partitions. This is the same
mathematics and supports arbitrary moderate dimension, though its cost grows
with the Bell number and is therefore not intended for very high dimension.

## Fitting API

The R package exposes many high-level fitting variants for fixed parameters,
nonstationary location design matrices, quantile reparameterizations and
profile-likelihood workflows. The Fortran API concentrates on stationary
numerical likelihood fitting and exposes lower-level likelihood/distribution
routines from which specialized constrained fits can be built. Those R
workflow variants are not duplicated as an object system.

Positive scales and bounded dependence parameters are optimized on smooth
unconstrained transforms, which is more robust with finite-difference BFGS
than optimizing invalid natural-parameter regions directly.

## Validation

`test_core.f90`, `test_bivariate_families.f90`, `test_reference.f90`, and `test_fit.f90` cover probability
inversions, independence identities, bivariate/multivariate consistency,
nonparametric bounds, POT likelihoods, clustering and optimization.

Reference values in `test_reference.f90` were evaluated directly from the
corresponding upstream C formulas for the same data and parameters.
