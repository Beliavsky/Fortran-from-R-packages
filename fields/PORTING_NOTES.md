# Porting notes

## Architecture

The translation is split into small numerical modules and an umbrella `fields` module. Dense linear algebra uses BLAS/LAPACK. Sparse Kriging uses the attached `spam-fortran` translation through its CSR matrix and sparse Cholesky APIs.

The upstream fixed-form file `fieldsF77Code.f` contains mature implementations of cubic smoothing splines, robust splines, polynomial basis construction, radial-basis multiplication, polygon membership, and Wendland-grid multiplication. Rather than replacing those algorithms, the source was mechanically converted to free form and placed in module `fields_native`; typed modern wrappers provide allocatable arrays, derived result types, explicit dimensions, and error checks.

## Kriging conventions

For covariance matrix `K`, smoothing/nugget ratio `lambda`, observation weights `w`, and drift matrix `T`, the dense fit uses

`M = K + lambda * diag(1/w)`.

It solves the universal-Kriging normal equations using SPD LAPACK routines, returns the random coefficients and drift coefficients, and computes the smoother trace, GCV, profile ML and profile REML quantities. Prediction covariance includes the uncertainty from estimating the drift. Sparse Kriging follows the same equations with spam Cholesky solves.

## Thin-plate splines

`Tps` is implemented using the radial basis plus polynomial null space and the augmented saddle system. `Tps.cov` follows the upstream cardinal-point construction. The exact smoother diagonal is retained because it is also needed by the translated `QTps` pseudo-data cross-validation calculation.

## Quantile splines

`QSreg` and `QTps` use the upstream iterative pseudo-response construction

`y_pseudo = f + C * psi_scale * psi((y-f)/psi_scale)`

with the asymmetric bounded score and loss from `qsreg.family.R`. If lambda is not supplied, the corresponding GCV smoothing-parameter search is repeated on each pseudo-response as in upstream `QSreg`.

## FFT and grids

The 2-D circulant embedding implementation uses a self-contained radix-2 complex FFT and increases the embedding until the covariance spectrum is nonnegative to numerical tolerance. `fft_interp_surface` mirrors upstream `interp.surface.FFT`; its direct spectral evaluation supports the odd source-grid dimensions required by that R routine without introducing an external FFT dependency for arbitrary target dimensions.

## Sparse dependency correction

The attached `spam-fortran` dependency had a wrapper defect in the `pivot='none'` branch of `spam_chol_factor`: `cholstepwise` uses `perm` immediately, but the wrapper had never initialized it. The vendored copy initializes `perm` and `invp` to the identity before the call. The MMD and RCM paths are unchanged.

## Native spline error flag

The inherited `css`/`rcss` routines only set `ierr` on failure and do not guarantee assignment of zero on success. Every modern wrapper therefore initializes `ierr=0` before calling these kernels. A nonlinear spline regression test was added because a linear-null-space-only test can accidentally hide this issue.

## Intentional API differences

Fortran routines accept arrays and explicit model parameters rather than R formulas, S3 objects, named lists, or `do.call` covariance functions. Where upstream permits arbitrary user R functions, this port exposes a set of built-in covariance families and lower-level covariance-matrix fitting APIs, so callers can construct a custom covariance matrix themselves.

`predictDerivative.Krig` is represented by `krig_predict_gradient_stationary`. Polynomial derivatives are analytic; the covariance component uses centered numerical differentiation so one implementation supports all translated stationary covariance families, including arbitrary-order Matérn.

The stochastic interchange/search presentation of `cover.design` is not reproduced; the deterministic minimax criterion and a greedy farthest-point design are supplied as the reusable numerical equivalents.
