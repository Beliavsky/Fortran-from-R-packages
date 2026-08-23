# Porting notes

## Upstream

This project translates the computational core of R package L1pack 0.62-4
(February 18, 2026), by Felipe Osorio and Tymoteusz Wolodzko. Upstream declares
`License: GPL-3`.

Original `DESCRIPTION`, `NAMESPACE`, R sources, and native C/Fortran sources are
retained under `orig/` for provenance.

## fastmatrix dependency

L1pack directly depends on fastmatrix. The previously translated
`fastmatrix-fortran` v0.2.0 tree is included under `vendor/` and referenced as an
FPM path dependency. Both projects are GPL-3.0-only, so the combined source tree
is license-compatible.

For the generalized spatial median, this project ports the original fastmatrix
AS 78 algorithm directly. This is intentional: the convenience `mediancenter`
routine in the earlier fastmatrix-fortran v0.2.0 translation computes
componentwise medians and is not a sufficient substitute for L1pack's use of the
true multivariate median center.

## Barrodale-Roberts L1 algorithm

`CALGO478.f` was translated from legacy fixed-form Fortran to free-form modern
Fortran with explicit declarations and module scoping. The simplex algorithm is
otherwise preserved.

As a parity check, the modern routine was run against the original upstream
`CALGO478.f` on the same nontrivial regression problem. Coefficients, residuals,
minimum absolute deviation, rank, status code, and iteration count matched to
floating-point precision.

## Multivariate Laplace Bessel weights

The Yavuz-Arslan EM update requires

`K_(p/2-1)(x) / K_(p/2)(x)`.

Because `p` is an integer, the needed Bessel orders are integer or half-integer.
The port uses scaled K0/K1 approximations plus recurrence for integer orders and
the exact finite half-integer expansion otherwise. This avoids an external
special-function dependency and avoids exponential underflow in the ratio.

## Multivariate fitting log likelihood

The upstream C fitter evaluates convergence using Mahalanobis distances computed
before the M-step while combining them with the newly updated scatter matrix.
The Fortran port recomputes the distances after each parameter update before
calculating the observed log likelihood. The returned `loglik` therefore equals
the sum of `log_dmlaplace` at the returned fitted parameters.

The same consistency correction is used by the generalized spatial-median fit.

## L1 concordance coefficient

The upstream Laplace CCC helper obtains a truncated first moment using R's
`integrate`. For its fixed Laplace scale, that integral has a closed form. The
Fortran implementation evaluates that exact expression directly; no numerical
quadrature is required.

Bootstrap variance estimation is retained and uses Fortran's intrinsic RNG.

## Plotting and R object machinery

Formula parsing, model frames, S3 dispatch, print/summary methods, diagnostic
plotting, and graphics-device code are omitted. Numerical counterparts such as
prediction, covariance/confidence intervals, simulation, quantile residuals,
and the simulated Laplace envelope are available directly through Fortran APIs.
