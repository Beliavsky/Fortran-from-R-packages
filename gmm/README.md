# gmm-fortran

A modern free-form Fortran/FPM translation of the computational core of the R package **gmm** 1.9-1 by Pierre Chausse.

The goal is numerical/statistical parity, not emulation of R's formula, S3, plotting, printing, or model-frame systems.

## Included computational functionality

- Generic nonlinear GMM from user-supplied moment and mean-Jacobian callbacks.
- Two-step, iterated, continuously updated (CUE), and fixed-weight GMM.
- MDS, iid, HAC, and identity moment covariance/weighting paths.
- Bartlett, Parzen, Tukey-Hanning, truncated, and quadratic-spectral kernels.
- Andrews automatic bandwidth and Wilhelm bandwidth calculations.
- Linear IV/GMM and 2SLS.
- System GMM plus SUR, 3SLS, FIVE-style, and random-effects numerical fits.
- Generalized empirical likelihood: EL, ET, CUE, HD, ETEL, ETHD, and RCUE.
- GEL Lagrange-multiplier solvers, implied probabilities, covariance calculations, and specification statistics.
- J tests and Kleibergen K/J/S weak-identification statistics.
- ATE GEL moment/Jacobian kernels for balance, sample-balance, and ATT formulations with linear, logit, and probit outcome models.
- ATE marginal effects and standard errors.
- Nolan parameterization 0/1 stable characteristic function (`char_stable`).
- Linear-algebra helpers required by the translated estimators.

`use gmm` imports the public numerical API.

## Minimal example

```fortran
use gmm
real(dp) :: y(8), x(8,2), z(8,3)
type(linear_gmm_result_t) :: fit
! ... fill y, x and z ...
call linear_gmm_fit(y,x,z,fit,method=LINEAR_TWO_STEP,covariance=COV_MDS)
print *, fit%coefficients
```

A complete example is in `example/linear_gmm_example.f90`.

## Build

With Fortran Package Manager:

```text
fpm build
fpm test
fpm run --example linear_gmm_example
```

The project uses standard free-form Fortran 2018. Linear solves, inverses, and
least-squares calculations use the local `rfortran-linalg` dependency and its
pinned pure-Fortran LAPACK backend. System BLAS and LAPACK are not required.

## Licensing

Code translated from **gmm** is licensed GPL-2.0-or-later, matching upstream `License: GPL (>= 2)`. The supplied `r_mod.f90` remains under the MIT License. See `LICENSE`, `LICENSES/`, `NOTICE.md`, and `upstream/`.

## Scope

See `API_MAPPING.md` for the correspondence with the R package and `PORTING_NOTES.md` for implementation choices and known interface differences.
