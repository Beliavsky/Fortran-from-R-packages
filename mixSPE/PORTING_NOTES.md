# Porting notes

## Upstream surface

The CRAN package exports `EMGr`, `rpe`, `rspe`, and an S3 print method. Internal
routines implement densities, covariance constraints, beta updates, skewness
updates, initialization and BIC model selection.

## Translation choices

1. R-only printing/S3 wrappers are omitted.
2. The internal elliptical PE density is exposed as `dpe`/`log_dpe`; upstream
   uses it internally and exposes PE simulation through `rpe`.
3. `dspe` is exposed because it is the numerical density underlying `rspe` and
   the skew-mixture E-step.
4. `mvtnorm-fortran` is vendored only for the source modules needed for
   multivariate-normal proposals and linear algebra.
5. The R package's covariance M-step contains specialized manifold optimization
   for several eigen-decomposed models. This port enforces the same covariance
   families by deterministic spectral projections of weighted scatter matrices.
   This is a generalized-EM implementation, not a statement of bitwise parity.
6. Semi-supervised labels are supported by `em_fit(labels=...)`, with zero for
   unknown and 1..G for known classes, as upstream.
