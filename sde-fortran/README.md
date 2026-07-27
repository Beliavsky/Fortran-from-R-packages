# sde-fortran

A modern Fortran translation and adaptation of the R package `sde` 2.0.21 by
Stefano Maria Iacus. The project provides simulation and inference routines for
one-dimensional stochastic differential equations (SDEs) in an FPM package.

The translation preserves the upstream GPL version 2 or later license. See
`LICENSE`, `NOTICE`, and the retained metadata in `upstream/`.

## Scope

The translated computational surface includes:

- Exact conditional and stationary laws for OU, geometric Brownian motion, and
  CIR models.
- Brownian motion, geometric Brownian motion, Brownian bridges, diffusion
  bridges, and exact-model path simulation.
- Euler predictor-corrector, Milstein, second-order Milstein, Kloeden-Platen-
  Schurz (KPS), Ozaki, Shoji, arbitrary conditional-law, and exact acceptance
  simulation.
- Euler, Elerian, Kessler, Ozaki, Shoji, Hermite, and Pedersen transition or
  likelihood approximations.
- Simple estimating functions, generator estimating functions, linear
  martingale estimating functions, two-stage GMM with Bartlett HAC weighting,
  Dacunha-Castelle-style likelihood/AIC, and divergence tests.
- Kernel drift, diffusion, and density estimators; diffusion changepoint
  detection; and Markov-operator distances based on B-splines.

R-specific expression parsing, `ts`/`zoo`/`fda` objects, plotting, console
formatting, and bundled R data files are not reproduced. Arbitrary model
expressions are represented by typed Fortran procedure callbacks.

See `COVERAGE.md` for a function-by-function map and `PORTING_NOTES.md` for
behavioral details and documented corrections.

## Requirements

- A Fortran 2018 compiler.
- Fortran Package Manager (FPM).

The project has no external library dependencies.

## Build and test

```text
fpm build
fpm test
fpm run sde_demo
```

Examples can be run with:

```text
fpm run --example custom_sde
fpm run --example model_laws
```

The tests exercise exact distribution round trips, all major simulation
families, density approximations, likelihoods, estimating functions, GMM,
AIC/divergence routines, nonparametric estimators, changepoints, B-splines,
and Markov-operator distances.

See `VALIDATION.md` for the compiler, strict flags, test results, and the
distinction between direct compiler validation and FPM invocation.

## Quick example

```fortran
program example_ou
   use sde
   implicit none

   real(dp), parameter :: theta(3) = [1.2_dp, 0.8_dp, 0.45_dp]
   real(dp), allocatable :: path(:, :)

   call seed_rng(12345_i64)
   call simulate_ou_exact([0.5_dp], 1.0_dp/252.0_dp, 500, theta, path)

   print '(a, f12.6)', "endpoint = ", path(size(path, 1), 1)
end program example_ou
```

The OU parameterization is

```text
dX(t) = (theta(1) - theta(2)*X(t)) dt + theta(3) dW(t).
```

The GBM parameterization is

```text
dX(t) = theta(1)*X(t) dt + theta(2)*X(t) dW(t).
```

The CIR parameterization is

```text
dX(t) = (theta(1) - theta(2)*X(t)) dt
        + theta(3)*sqrt(X(t)) dW(t).
```

## Custom coefficients

General simulators use callbacks with this interface:

```fortran
pure function coefficient(t, x, theta) result(value)
   use sde_kinds, only : dp
   real(dp), intent(in) :: t, x, theta(:)
   real(dp) :: value
end function coefficient
```

For example:

```fortran
call simulate_euler([0.4_dp], 0.0_dp, 0.002_dp, 500, drift, diffusion, &
   theta, path, predictor_corrector=.false.)
```

Time-independent likelihood and information routines use:

```fortran
pure function state_function(x, theta) result(value)
   use sde_kinds, only : dp
   real(dp), intent(in) :: x, theta(:)
   real(dp) :: value
end function state_function
```

Additional callback interfaces for transition samplers, moment functions,
estimating functions, generator derivatives, and martingale weights are in
`src/sde_interfaces.f90`.

## Path layout

Multi-path simulators return an allocatable array with shape
`(n_steps + 1, n_paths)`. The first row contains the supplied initial values.
Single-path Brownian and bridge routines return rank-one arrays of length
`n_steps + 1`.

## Modules

- `sde`: umbrella module for the public API.
- `sde_models`: exact OU, GBM, and CIR laws.
- `sde_simulation`: path simulators and bridges.
- `sde_density`: local transition-density approximations.
- `sde_likelihood`: Euler, Hermite, Pedersen, and exact-model likelihoods.
- `sde_estimating`: estimating-function and martingale methods.
- `sde_gmm`: iterative two-stage GMM with HAC weighting.
- `sde_information`: Dacunha-Castelle-style log likelihood and AIC.
- `sde_divergence`: likelihood-ratio and phi-divergence tests.
- `sde_nonparametric`: kernel drift, diffusion, and density estimators.
- `sde_change_point`: diffusion-scale changepoint estimation.
- `sde_markov_distance`: B-spline Markov-operator distances.
- `sde_special`, `sde_distributions`, `sde_random`, `sde_linalg`, and
  `sde_optimization`: dependency-free numerical support.

## Numerical design

The library uses `dp = kind(1.0d0)`, `implicit none`, explicit interfaces,
allocatable results, procedure callbacks, and no global model state. Random
number generation uses the compiler intrinsic generator with local normal,
gamma, Poisson, chi-square, and noncentral chi-square transforms. Matrix
operations and bounded Nelder-Mead optimization are implemented internally.

For production inference, inspect convergence flags and numerical diagnostics,
choose discretization and Monte Carlo controls appropriate to the sampling
interval, and validate results against model-specific theory or an independent
implementation.

## License

`sde-fortran` is a derivative work of `sde` 2.0.21 and is distributed under
GNU GPL version 2 or any later version (`GPL-2.0-or-later`). Every translated
source file carries an SPDX identifier. The complete GPL version 2 text is in
`LICENSE`.
