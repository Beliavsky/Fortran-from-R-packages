# fracdiff-fortran

A modern Fortran 2018 and FPM translation of the computational algorithms in
R package `fracdiff` 1.5-4.

The project implements the complete numerical scope of the package, not only
the `fracdiff()` functionality required by `ufRisk`.

## Implemented algorithms

- Haslett-Raftery approximate Gaussian maximum likelihood for ARFIMA(p,d,q)
- Outer Brent minimization over the fractional parameter `d`
- Inner Levenberg-Marquardt nonlinear least-squares estimation of AR and MA terms
- Analytical ARMA residual Jacobian
- Numerical likelihood Hessian, covariance, standard errors, and correlations
- Covariance recomputation with a user-specified finite-difference interval
- Exact fractional-noise simulation followed by inverse ARMA filtering
- User-supplied innovation and burn-in sequences
- Fast FFT fractional differencing and the original direct convolution method
- Geweke-Porter-Hudak estimator
- Smoothed-periodogram/Sperio estimator
- Coefficient extraction, confidence intervals, AIC, BIC, and summary statistics
- Polynomial multiplication and AR-polynomial root calculations used by simulation

The MA convention is the one used by the original package:

```text
x(t) - ar(1)*x(t-1) - ... = e(t) - ma(1)*e(t-1) - ...
```

Thus the signs of the MA coefficients are reversed relative to the convention
used by many other ARMA implementations.

## Build with FPM

```sh
fpm build
fpm test
fpm run
fpm run --example long_memory_estimators
fpm run --example custom_innovations
```

The project has no external library dependencies.

## Minimal example

```fortran
program example
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff, only : dp, fracdiff_simulation, fracdiff_model, &
                        fracdiff_sim, fracdiff_fit
   implicit none

   type(fracdiff_simulation) :: simulated
   type(fracdiff_model) :: fitted

   simulated = fracdiff_sim(2000, 0.3_dp, ar=[0.2_dp], ma=[-0.4_dp], &
                            n_start=100, seed=107_int64)
   fitted = fracdiff_fit(simulated%series, nar=1, nma=1)

   print *, fitted%d
   print *, fitted%ar
   print *, fitted%ma
   print *, fitted%sigma
end program example
```

## Project layout

- `src/`: Fortran library
- `app/`: command-line demonstration
- `example/`: focused examples
- `test/`: numerical and structural tests
- `original/`: original R package source retained for provenance
- `API.md`: public API reference
- `PORTING.md`: detailed source-to-Fortran mapping and intentional differences
- `TESTING.md`: validation procedure and coverage
- `REFERENCE_GENERATION.md`: independent reference calculations

## Licensing

The original package declares `GPL (>= 2)`. This translation is therefore
licensed under `GPL-2.0-or-later`. Complete GPL-2.0 and GPL-3.0 texts are in
`LICENSES/`.
