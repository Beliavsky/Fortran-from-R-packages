# Porting notes

## Scope

This project translates the computational/statistical layer of R `urca`
1.3-4 to modern Fortran. It does not attempt to reproduce R's S4 classes,
formula evaluation, plotting, printing, or dataset loading.

## Linear algebra

The R implementation relies heavily on `lm`, `solve`, `chol`, `eigen`, and
QR calculations. The Fortran translation uses BLAS/LAPACK through small typed
wrappers. Core calls include linear least squares, LU/SPD inversion,
Cholesky factorization, symmetric eigensystems, and QR construction.

Johansen estimation is performed through residual covariance matrices and a
symmetric whitened eigenproblem. Eigenvectors are normalized in the same
logical convention used by the R routines before loading/PI calculations.

## Johansen restrictions

`blrtest`, `alrtest`, `ablrtest`, `bh5lrtest`, and `bh6lrtest` were translated
from the upstream matrix algorithms. In particular, `bh6lrtest` preserves an
upstream-specific iterative update in which the raw constrained eigenvectors
are used when rebuilding the known subspace. Replacing that step with the
whitening matrix changes the statistic measurably.

## Zivot-Andrews indexing

The lagged-difference block is aligned to reproduce the R `embed`/data-frame
construction. During validation an initial one-observation shift in this
translation was detected and corrected. The deterministic reference test now
reproduces the independently reconstructed statistic and break index.

## Phillips-Perron auxiliary statistics

The port returns the auxiliary statistics exposed by upstream `ur.pp`, not
only the headline Z-alpha or Z-tau statistic. They are stored in
`ur_test_result%auxiliary_statistics`.

## MacKinnon response surfaces

The upstream package contains a legacy fixed-form Fortran implementation and
large response-surface tables represented in R. The port:

1. retains the original files unchanged under `upstream/`;
2. converts the response-surface data to Fortran parameter arrays;
3. reimplements the small GLS/evaluation routines in modern Fortran;
4. keeps the parameter arrays in groups small enough to remain within
   standard Fortran continuation limits.

The upstream package's MacKinnon license/permission correspondence remains at
`upstream/inst/Licenses/MacKinnonLicense.txt`.

## `nlme` dependency audit

Although upstream NAMESPACE contains `import(nlme)`, a source audit found no
`nlme` routine used by the exported numerical algorithms translated here.
The user-supplied `nlme-fortran` project is therefore retained under
`reference/nlme-fortran/` as a compatibility/reference artifact but is not a
link dependency.

## Numerical status codes

Fortran routines expose integer `info` fields/arguments instead of throwing R
conditions. Zero indicates normal completion unless a routine documents an
iterative nonzero informational status. Negative values generally represent
invalid dimensions/options; positive values generally propagate numerical
factorization/solve failures or iteration status.

## Source format

All translated sources are free-form Fortran 2018, use `implicit none` (via
module/project policy), explicit `real(dp)` precision, allocatable arrays,
and explicit LAPACK interfaces. No fixed-form source is compiled by the FPM
project.
