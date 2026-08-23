# Porting notes

## Numerical special functions

The R package obtains arbitrary-order modified Bessel K and incomplete Bessel K
functionality from R/DistributionUtils. This standalone port evaluates K_nu(x)
from its integral representation, with a large-x asymptotic branch. GIG and GH
CDFs are evaluated by numerical quadrature and quantiles by safeguarded
bisection. This is portable and dependency-free, but slower than the mature R
special-function implementations, especially for repeated tail probabilities.

## Fitting

`gig_fit`, `hyperb_fit`, and `nig_fit` preserve the corresponding likelihoods
and parameter constraints but use deterministic coordinate-pattern searches
instead of R's `nlm`, `optim`, and `nlminb`. Therefore fitted values should be
numerically comparable but are not expected to reproduce R optimizer traces or
last-bit parameter values.

`hyperblm_fit` provides the numerical linear-model entry point and hyperbolic
residual fit, but omits R formula/model-frame machinery and S3 diagnostics.

## Deliberately omitted R-facing functionality

Plot/QQ helpers, S3 print/summary/vcov wrappers, datasets, and graphics are not
part of the Fortran library. The Cramer-von Mises presentation/test wrapper is
also omitted from v0.1.0; its distribution primitives are available for a
future dedicated inference module.
