# mcmc-fortran

Modern Fortran/FPM translation of the computational core of the R package
`mcmc` 0.9-8 by Charles J. Geyer and Leif T. Johnson.

Implemented:

- random-walk Metropolis MCMC (`metrop`);
- scalar, diagonal, and full-matrix proposal scales;
- batch means and per-batch acceptance rates;
- optional output callbacks and debug trajectories;
- serial simulated tempering / umbrella sampling;
- parallel tempering with swaps;
- morphometric transformations and `morph_metrop`;
- initial-positive, initial-monotone, and initial-convex sequence variance
  estimators (`initseq`);
- overlapping batch means covariance (`olbm`).

The library uses typed Fortran callbacks rather than R functions/environments
and has no runtime dependencies.

## Minimal Metropolis example

```fortran
use mcmc

type(metrop_result) :: fit

call set_mcmc_seed(1234)
fit = metrop(log_target, [0.0_dp,0.0_dp], nbatch=500, blen=20, &
             scale=scale_constant(1.0_dp))
```

A callback has the form

```fortran
subroutine log_target(state,value,data)
   use mcmc, only : dp
   real(dp), intent(in) :: state(:)
   real(dp), intent(out) :: value
   class(*), intent(in), optional :: data
   value = -0.5_dp*sum(state*state)
end subroutine
```

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The complete supplied R/C source is retained under `upstream/`.
