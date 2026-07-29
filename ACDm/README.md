# ACDm for Modern Fortran

A self-contained modern Fortran/FPM translation of the computational algorithms
in the R package **ACDm 1.1.0**.

ACDm provides estimation, simulation, distribution support, transaction-duration
construction, diurnal adjustment, residual diagnostics, and specification tests
for autoregressive conditional duration models.

## Numerical scope

The library implements all 13 duration dynamics from ACDm:

- ACD
- LACD1 and LACD2
- EXACD
- AMACD
- ABACD and AACD
- TACD and TAMACD
- BACD and BCACD
- SNIACD and LSNIACD

It also implements the package's ten duration distributions:

- exponential
- Weibull
- Burr
- generalized gamma
- generalized F
- q-Weibull
- q-Weibull/exponential mixture
- q-Weibull/Weibull mixture
- finite inverse-Gaussian mixture
- Birnbaum-Saunders

Major facilities include:

- maximum-likelihood estimation with fixed parameters and box bounds
- numerical Hessian, score matrix, conventional and robust covariance estimates
- parametric bootstrap standard errors
- deterministic conditional forecasts
- simulation with generated or user-supplied innovations
- trade, price, and volume duration construction
- cubic-spline, penalized-spline, adaptive super-smoother, and Flexible Fourier
  Form diurnal adjustment
- probability-integral and Cox-Snell residual transforms
- ACF, QQ, density, hazard, and likelihood-profile numerical diagnostics
- Meitz-Terasvirta remaining-ACD, STACD, and TVACD LM tests

The implementation has no mandatory external numerical dependency.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example diurnal_adjustment
fpm run --example threshold_model
```

A reproducible GNU Fortran script is also included:

```text
./run_gfortran_tests.sh strict
./run_gfortran_tests.sh optimized
```

On Windows with GNU Fortran:

```text
run_gfortran_tests.bat strict
run_gfortran_tests.bat optimized
```

## Minimal example

```fortran
program example
  use acdm
  implicit none

  integer, parameter :: n = 500
  type(acd_order) :: order
  type(acd_fit_options) :: options
  type(acd_fit_result) :: fit
  type(rng_state) :: rng
  real(dp) :: durations(n)
  integer :: status

  order = acd_order(p=1, r=0, q=1)
  call seed_rng(rng, 12345)
  call simulate_acd(n, MODEL_ACD, order, [0.2_dp, 0.15_dp, 0.7_dp], &
                    DIST_EXPONENTIAL, [real(dp) ::], durations, status, rng)

  options%model = MODEL_ACD
  options%dist = DIST_EXPONENTIAL
  options%order = order
  call acd_fit_model(durations, options, fit)

  print *, fit%parameters
end program example
```

## Modules

- `acdm_kinds`: kinds and status constants
- `acdm_math`: RNG, special functions, quantiles, linear algebra, splines
- `acdm_distributions`: duration distributions and generic distribution API
- `acdm_models`: model recursions, likelihoods, and simulation
- `acdm_fit`: estimation, inference, bootstrap, and forecasting
- `acdm_data`: duration construction and diurnal adjustment
- `acdm_diagnostics`: residual diagnostics and LM tests
- `acdm_profiles`: hazard estimates and likelihood profiles
- `acdm`: primary aggregate module
- `acdm_api`: aliases resembling the original R entry-point names

## Intentional omissions

The following R-runtime or presentation features are not numerical algorithms and
are not reproduced:

- S3 classes and methods
- `broom` tidiers
- ggplot2/rgl graphics
- interactive trace plots
- R data-frame and POSIXlt dispatch
- bundled serialized demonstration datasets
- R optimizer adapters

The numeric information used by the main diagnostic plots is available through
Fortran result types.

## Documentation

- `API.md`: procedure and result-type reference
- `PORTING.md`: source mapping and compatibility decisions
- `TESTING.md`: validation and build details
- `REFERENCE_GENERATION.md`: independent-reference methodology
- `CHANGELOG.md`: changes introduced by the translation

## License

GPL-3.0-or-later. See `LICENSE`, `NOTICE.md`, and the retained original package.
