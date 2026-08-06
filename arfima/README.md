# arfima-fortran

A modern Fortran 2018 translation of the computational code in the R package
`arfima` 1.8-2 by J.Q. Veenstra and A.I. McLeod.

The library fits, simulates, filters, and forecasts long-memory time series
that combine ordinary and seasonal ARMA terms with one of:

- fractionally differenced white noise (FDWN),
- fractional Gaussian noise (FGN),
- power-law autocovariance (PLA) noise, or
- no long-memory component.

## Main features

- Exact ARMA, FDWN, FGN, PLA, FARMA, and seasonal ARFIMA autocovariances
- Exact Gaussian likelihood and innovations through Durbin-Levinson recursion
- Stable AR/MA parameterization through partial autocorrelations
- Integer and seasonal differencing and inverse differencing
- Static regressors and Box-Jenkins dynamic transfer functions
- Bounded long-memory parameter transformations
- Maximum-likelihood fitting and multiple-start mode searches
- Numerical Hessian covariance estimates, AIC, and BIC
- Exact conditional Gaussian forecasts, including integrated models
- Simulation from theoretical autocovariances
- Numerical Fisher information and identifiability checks
- R-name compatibility wrappers for the main computational entry points

Plotting, R formula/S3 infrastructure, parallel clusters, package datasets, and
formatted print methods are omitted. The original package archive is retained
under `upstream/`.

## Build

With GNU Make and GNU Fortran:

```sh
make check
make optimized
make example
```

With FPM:

```sh
fpm test
fpm run --example arfima_demo
```

The checked build uses bounds, allocation, floating-point, and uninitialized
value diagnostics. The optimized build uses `-O3` with warnings treated as
errors.

## Minimal example

```fortran
program demo
  use arfima
  implicit none

  type(arfima_spec) :: spec
  type(arfima_parameters) :: truth
  type(arfima_fit_result) :: fit
  type(arfima_error) :: error
  real(dp), allocatable :: x(:)

  spec%p = 1
  spec%lmodel = long_memory_fd

  allocate(truth%phi(1), truth%theta(0), truth%phiseas(0), truth%thetaseas(0))
  allocate(truth%beta(0), truth%delta(0), truth%omega(0))
  truth%phi = 0.4_dp
  truth%dfrac = 0.2_dp
  truth%mean = 0.5_dp

  call set_random_seed(12345)
  call arfima_simulate(spec, truth, 200, 1.0_dp, x, error)
  call fit_arfima(spec, x, fit)

  write(*,*) fit%parameters%phi
  write(*,*) fit%parameters%dfrac
end program demo
```

## Numerical conventions

Moving-average coefficients use the Box-Jenkins sign convention used by the R
package: the MA polynomial is `1 - theta(1) B - ...`.

The likelihood is the concentrated exact Gaussian likelihood up to the same
additive constant omitted by the upstream package. `sigma2` is estimated from
the Durbin-Levinson innovations.

## License

MIT, matching the upstream package. See `LICENSE` and
`licenses/UPSTREAM_LICENSE`.
