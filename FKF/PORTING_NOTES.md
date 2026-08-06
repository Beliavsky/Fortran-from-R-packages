# Porting notes

## Numerical formulation

The filter follows the upstream covariance-form recursion. Positive-definite
innovation covariance matrices are factored by a self-contained lower-Cholesky
routine, from which both the inverse and log determinant are obtained. State
covariances are symmetrized after updates to remove roundoff asymmetry.

The smoother implements the same Durbin-Koopman disturbance recursion as the C
source:

- `L_t = T_t - T_t K_t Z_t`
- `r_{t-1} = Z_t' F_t^-1 v_t + L_t' r_t`
- `N_{t-1} = Z_t' F_t^-1 Z_t + L_t' N_t L_t`

## Missing-data likelihood constant

The upstream C routine initializes the Gaussian constant with `n*d` even when
some or all components of an observation are missing. Consequently, each
missing scalar still contributes `-0.5*log(2*pi)` to the returned likelihood.

For source parity this behavior is the default. Pass
`corrected_missing_likelihood=.true.` to count only observed scalar values.
The state estimates are identical in both modes.

## Interface changes

R closures, S3 classes, plotting, elapsed-time fields, and R memory management
are not reproduced. Typed Fortran model/result objects and explicit status
codes are used instead. The implementation has no external BLAS/LAPACK runtime
dependency.
