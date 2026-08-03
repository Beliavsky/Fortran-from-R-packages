# Porting notes

## Data orientation

The R convention is preserved: rows are observations and columns are variables.

## Missing observations

R uses `NA`. The Fortran API uses IEEE quiet NaNs. Tyler, Cauchy, and skewed-t
fits drop incomplete rows. `fit_mvt(...,na_rm=.false.)` retains them and applies
the conditional first- and second-moment E-step from the original package.

## Linear algebra

The translation is dependency-free. It includes dense pivoted solves,
Cholesky log determinants, a Jacobi symmetric eigensolver, and the factor-model
update required by `fit_mvt`. These are intended for small and medium dense
problems, matching the original package's typical use. Large applications may
benefit from replacing the internal kernels with BLAS/LAPACK calls.

## Special functions

Fortran supplies `log_gamma` but not portable real-order modified Bessel K or
digamma functions. The translation therefore includes:

- an asymptotic-recurrence digamma implementation;
- scaled adaptive quadrature for fractional-order Bessel K base values;
- log-domain order recurrence for large real orders;
- finite-difference order derivatives for the skewed-t E-step.

This avoids a mandatory external special-function library.

## Optimization

R's `optimize()` and `uniroot()` calls are replaced by bounded golden-section
search and bracketed scalar root iteration. The resampled diagonal MLE and the
simulation-based POP correction use a deterministic portable random-number
generator so tests are reproducible.

## Result representation

R lists become the `heavy_tail_fit` derived type. Iteration-history lists used
only for plotting are not stored. Final parameters, convergence information,
CPU time, log likelihood, latent weights, and factor-model components remain
available.

## Gaussian limit

As in the R package, very large `nu` values are used as a numerically stable
Gaussian limit. `fixed_nu=1.0e12_dp` reproduces the maximum-likelihood sample
mean and covariance to numerical precision.

## Naming

R dots and capitalization are converted to lower-case underscore names, such
as `fit_Tyler` to `fit_tyler` and `nu_POP_estimator` to
`nu_pop_estimator`. Method-selection strings retain their R spelling.
