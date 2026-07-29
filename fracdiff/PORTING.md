# Porting notes

## Functional mapping

| R package entry point | Fortran API |
|---|---|
| `fracdiff()` | `fracdiff_fit()` |
| `fracdiff.var()` | `fracdiff_var()` |
| `fracdiff.sim()` | `fracdiff_sim()` |
| internal `fdsim` | `fractional_arma_filter_simulation()` |
| `diffseries()` | `diffseries()` |
| internal `diffseries0()` | `diffseries_direct()` |
| `fdGPH()` | `fd_gph()` |
| `fdSperio()` | `fd_sperio()` |
| `coef.fracdiff()` | `fracdiff_coefficients()` |
| `vcov.fracdiff()` | `model%covariance` |
| `residuals.fracdiff()` | `model%residuals` |
| `fitted.fracdiff()` | `model%fitted` |
| `logLik.fracdiff()` | `model%log_likelihood` |
| `confint.fracdiff()` | `fracdiff_confint()` |
| `summary.fracdiff()` | `summarize_fracdiff()` |
| internal `poly_mult` | `polynomial_multiply()` |

## Maximum-likelihood algorithm

The port retains the package's two-level estimation design:

1. The Haslett-Raftery fractional filter computes the approximate innovation
   sequence, innovation-variance terms, and generalized mean estimate.
2. For a fixed `d`, AR and MA coefficients minimize the conditional residual
   sum of squares by Levenberg-Marquardt iteration using the source analytical
   residual Jacobian.
3. Brent minimization profiles the likelihood over `d`.
4. The Gaussian profile log likelihood uses the source constant `2.8378` and
   the same innovation-variance divisor.

## Numerical-library modernization

The original implementation was old Fortran translated to C with `f2c`, plus
MINPACK, LINPACK, SLATEC gamma routines, and R's BLAS. The port is standalone
modern Fortran:

- a compact Levenberg-Marquardt implementation replaces the bundled MINPACK
  translation while retaining the same residual vector and analytical Jacobian;
- intrinsic `gamma` replaces the SLATEC gamma implementation;
- partial-pivoting linear solves replace LINPACK SVD inversion;
- the covariance matrix is based on a symmetric finite-difference Hessian of
  the complete fixed-parameter likelihood;
- a radix-2 FFT replaces R's `fft()` and `nextn()` calls.

These substitutions remove external dependencies and global common blocks.
They can lead to small numerical differences in optimizer paths, covariance
estimates, and failure handling, but do not remove any package-level algorithm.

## Simulation

The fractional-noise recursion and inverse ARMA filter follow `fdsim` directly.
The source burn-in and `backComp` slicing rules are retained. Differences are:

- the default random generator is portable Park-Miller plus Box-Muller rather
  than R's RNG stream;
- callers can supply `innovations` and `start_innovations` for exact,
  deterministic comparisons;
- the mean is added after ARMA filtering, as in the original low-level routine.

## Residual and fitted-value initialization

The R wrapper calls `stats::arima()` after estimation solely to construct the
reported residual and fitted series. The Fortran port uses explicit conditional
ARMA recursion with unavailable pre-sample values set to zero. This does not
alter the Haslett-Raftery fit or likelihood, but the earliest reported residuals
may differ from R's state-space/Kalman initialization.

## R-only infrastructure omitted

The following are presentation or runtime features rather than numerical
algorithms and are not reproduced:

- S3 class registration and dispatch
- formatted `print()` methods
- R `ts` attributes and time bases
- `NA` handling conventions
- arbitrary R `rand.gen` callbacks
- R environment-variable test controls

Typed result objects, explicit status codes, supplied innovations, and normal
Fortran I/O provide the corresponding functionality.
