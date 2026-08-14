# Algorithm notes

## Numerical architecture

`ks-fortran` evaluates Gaussian kernels directly at arbitrary points and uses a
single derivative-tensor implementation for KDDE and derivative functionals.
Grid acceleration is separated into generic linear binning/interpolation and
symmetric convolution routines rather than maintaining four dimension-specific
C entry points.

SPD matrix operations use LAPACK through explicit Fortran interfaces.  This is
intentional: FPM is configured with `implicit-external = false`, so accidental
implicit BLAS/LAPACK interfaces are compile-time errors.

Random generation uses an explicit local RNG state instead of R's global RNG.
This makes multiple independent streams possible and avoids hidden global state.

## Bandwidth selection

The 1D normal-reference, LSCV, staged plug-in and SCV paths follow the ks
criteria.  The multivariate `hpi_matrix` in v0.1.0 is a self-contained one-stage
plug-in implementation using a fourth-derivative pilot functional and bounded
positive-definite parameterization.  It is not a claim that every R `Hpi`
pilot/preconditioning/optimizer option has been reproduced.  The full `Hscv`,
`Hbcv`, and normal-mixture (`Hnm`) selector families are retained as explicit
future-work items rather than silently aliased to another selector.

KDA defaults to a normal-reference bandwidth if `Hs` is not supplied; the R
package's higher-level `Hkda` wrappers can select more elaborate plug-in
bandwidths.  Callers needing exact selector control should supply `Hs`.

## KCDE and mvtnorm

`kcde` is a kernel **cumulative distribution** estimator.  For dimension > 1,
each Gaussian component CDF is evaluated with the vendored `mvtnorm-fortran`
Genz/Bretz probability implementation.  The supplied mvtnorm port is
GPL-2.0-only, so this combined distribution elects the GPL-2 option offered by
upstream ks.

The upstream exact-point R helper accepts a weight argument but does not apply
it in its `pmvnorm` loop.  The Fortran `kcde_eval` applies normalized weights;
this is an intentional correctness fix and is consistent with weighted KDE
semantics elsewhere in ks.

## Deconvolution

The R package uses `kernlab::ipop` for a simplex-constrained quadratic program.
The Fortran port solves the same quadratic objective with a native projected-
gradient simplex solver, keeping the library independent of an external QP
package.

## Density ridges

`density_ridge_point` follows the Gaussian mean-shift form of upstream
`kdr.base`: the mean-shift vector is projected onto the Hessian eigenvectors
normal to the requested ridge and iterated.  It does not replace that update by
an arbitrary gradient step.

## Intentional interface omissions

R formula/S3 dispatch, plotting, colors, progress bars, data-frame validation,
and display-only contour/perspective machinery are omitted.  Computational
branches listed as deferred in `API_MAPPING.md` are also not represented by
misleading stubs.
