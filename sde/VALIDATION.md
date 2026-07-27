# Validation

## Environment

The supplied project was validated on 2026-07-25 with:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

The FPM executable was not installed in the translation environment. The
`fpm.toml` file was parsed as TOML and checked against the current FPM manifest
layout. The complete source graph, tests, application, and examples were then
compiled and linked directly with GNU Fortran in module-dependency order.

## Strict compiler flags

```text
-std=f2018
-O0
-g
-Wall
-Wextra
-Wimplicit-interface
-Wconversion-extra
-Werror
-fcheck=all
-fbacktrace
-ffree-line-length-none
```

## Results

All four test programs passed:

```text
test_models: PASS
test_simulation_density: PASS
test_inference: PASS
test_nonparametric: PASS
```

The application and both examples also compiled, linked, and ran:

```text
sde_demo
custom_sde
model_laws
```

GNU ld reports executable-stack warnings for objects that pass internal
Fortran procedures as callbacks. This is a GNU trampoline implementation
detail; it did not produce a compiler warning or test failure. Applications
with a strict non-executable-stack policy can replace internal callbacks with
module procedures or adapt the callback/context interfaces.

## FPM commands for users and CI

From the project root:

```text
fpm build
fpm test
fpm run sde_demo
fpm run --example custom_sde
fpm run --example model_laws
```

The default `src`, `app`, `example`, and `test` directories are used, with
automatic target discovery enabled in `fpm.toml`.

## Numerical scope of validation

The tests cover:

- OU, GBM, and CIR conditional and stationary distribution identities.
- Random generators, CDF/quantile round trips, and exact-model likelihoods.
- Brownian, bridge, exact-model, Euler, Milstein, KPS, Ozaki, and Shoji paths.
- Euler, Elerian, Kessler, Ozaki, Shoji, Hermite, and Pedersen approximations.
- Estimating functions, GMM, numerical optimization, AIC, and divergence tests.
- Kernel estimators, changepoints, B-spline bases, missing-value interpolation,
  and Markov-operator distances.

Monte Carlo checks use fixed seeds and statistical tolerances. They establish
basic implementation consistency, not bit-for-bit equivalence with R or a
complete numerical certification over every parameter regime.
