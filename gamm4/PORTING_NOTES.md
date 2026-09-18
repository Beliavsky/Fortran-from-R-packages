# Porting notes

## Why this package is a bridge

There is little stand-alone statistical machinery in upstream `gamm4`. Its R
code asks mgcv to construct smooths, converts penalized smooth coefficients to
random effects, calls lme4, then reconstructs a GAM-like covariance/EDF view.
The Fortran package preserves that architecture through sibling FPM
dependencies.

## Smooth reparameterization

For a single quadratic penalty `S`, the translation computes a symmetric
eigendecomposition. Eigenvectors in the numerical null space become unpenalized
fixed-effect columns. Positive-eigenvalue vectors are divided by the square
root of their eigenvalues, giving an identity penalty and hence one isotropic
random-effect variance per smooth. The inverse relation is the smoothing
parameter `sp = 1 / tau^2` on the relative mixed-model scale.

Multi-penalty smooths are rejected rather than silently collapsed to one
penalty. Supporting those terms faithfully requires the more general
`smooth2random` block bookkeeping used by mgcv/gamm4.

## Mixed-model objective

The Gaussian path profiles the residual scale and implements ML/REML criteria
using the same relative-covariance interpretation as lme4. The generalized
path uses penalized IRLS conditional modes and the lme4-style Laplace
log-determinant correction.

Ordinary random-effect covariance matrices are created through the translated
lme4 covariance parameterizations, so grouped intercept/slope terms can use
unstructured, diagonal, compound-symmetry, or AR(1) structures.

## Coefficient covariance (`getVb`)

The direct upstream formula is retained:

`V = diag(v) + scale * Z * phi * Z'`

and the penalized fitting-space covariance is transformed back to the original
GAM coefficient parameterization. The upstream Python/CHOLMOD route is an
acceleration of this calculation and is therefore omitted.

## Interface differences

R formulas, model frames, factors, contrasts, S3 class construction, and
`reticulate` are intentionally absent. Numeric arrays make the statistical
operation explicit and keep the Fortran API independent of an R runtime.
