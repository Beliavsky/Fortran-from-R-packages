# rumidas-fortran

Modern Fortran translation of the computational core of the R package
`rumidas` 0.1.3.

The library implements univariate GARCH-MIDAS, double-asymmetric
GARCH-MIDAS, GARCH-MIDAS-X, two-MIDAS-variable models, MEM, MEM-X,
MEM-MIDAS, and MEM-MIDAS-X.  It uses the attached modern Fortran
translation of `maxLik` as an FPM path dependency.

## Build

```text
fpm build
fpm test
fpm run
```

GNU Fortran validation without FPM is available through:

```text
sh scripts/test_gfortran.sh
sh scripts/test_gfortran_optimized.sh
```

On Windows with `gfortran` in `PATH`:

```text
scripts\test_gfortran.bat
```

## Public module

```fortran
use rumidas
```

All real calculations use:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

## Basic GARCH-MIDAS evaluation

```fortran
use rumidas
implicit none

type(garch_midas_spec) :: spec
real(dp) :: parameters(5)
real(dp) :: returns(100), midas_matrix(13, 100)
real(dp), allocatable :: loglik(:), volatility(:), long_run(:), short_run(:)
integer :: status

spec = garch_midas_spec(RUMIDAS_GM, RUMIDAS_NORMAL, &
  RUMIDAS_BETA_LAG, 12, 0, .false.)
parameters = [0.05_dp, 0.90_dp, -8.0_dp, 0.20_dp, 3.0_dp]

call garch_midas_evaluate(parameters, spec, returns, midas_matrix, &
  loglik, volatility, long_run, short_run, status)
```

The MIDAS matrix has shape `(K+1, number_observations)`.  Its columns contain
the low-frequency observations from lag `K` through the contemporaneous
period.  As in the R implementation, the contemporaneous row receives zero
weight, preventing look-ahead.

## Estimation

```fortran
type(rumidas_fit_result) :: fit
type(rumidas_fit_control) :: control

control%random_starts = 20
control%method = 'bfgs'
call ugmfit(spec, returns, midas_matrix, fit, status, control=control)
```

`fit` contains coefficients, Hessian and robust standard errors, fitted total,
long-run, and short-run variances, observation likelihoods, AIC, BIC, and
optimizer diagnostics.

## Scope

Direct numerical coverage includes:

- normalized Beta and exponential-Almon MIDAS weights;
- explicit period-index construction of MIDAS lag matrices;
- GM, GMX, GM2M, DAGM, DAGMX, and DAGM2M;
- skewed and non-skewed short-run recursions;
- Gaussian and variance-standardized Student-t likelihoods;
- MEM, MEMX, MEMMIDAS, and MEMMIDASX;
- all upstream likelihood, conditional prediction, and long-run helper paths;
- deterministic random-start selection and `maxLik` estimation;
- Hessian and observation-score sandwich covariance estimates;
- MSE, QLIKE, AIC, BIC, and multi-step variance forecasts.

R `xts`/`zoo` dates, S3 printing, packaged data, and automatic time-frequency
aggregation are intentionally not reproduced.  See `API_MAP.md` and
`PORTING_NOTES.md`.

## Licensing

The original `rumidas` package is GPL-3.  This translation is distributed under
`GPL-3.0-only`.  The vendored `maxLik` translation is GPL-2.0-or-later and is
compatible with this combined distribution.  The complete upstream source is
retained under `original/rumidas-master`.
