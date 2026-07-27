# Validation

## Toolchain

Validation was performed with GNU Fortran 14.2.0 in Fortran 2018 mode.

The library, tests, application, and examples were compiled with:

```text
-std=f2018
-Wall
-Wextra
-Wconversion-extra
-Wimplicit-interface
-Werror
-fcheck=all
-fbacktrace
```

## Test programs

The following programs compiled and passed:

```text
test_core: PASS
test_bootstrap: PASS
test_spectral: PASS
```

`test_core` checks known normal and chi-square quantiles, deterministic AR(1),
variance-ratio, Lo-MacKinlay, Chow-Denning, Wald, Wright, and curve outputs.
The deterministic AR(1), VR(2), M1, and M2 values were independently reproduced
with NumPy formulas to machine precision.

`test_bootstrap` exercises wild and automatic bootstraps, panel VR,
subsampling, and single/joint Wright critical-value simulation. It checks
array dimensions, probability bounds, and ordered empirical quantiles.

`test_spectral` exercises the automatic portmanteau, average exponential,
spectral shape, generalized spectral, Dominguez-Lobato, and Chen-Deo tests,
including the Chen-Deo covariance solve.

## Executable smoke tests

These targets compiled and ran successfully:

```text
vrtest_demo
basic_tests
bootstrap_tests
```

## FPM status

The runtime used for translation did not have the `fpm` executable installed,
so direct `fpm build` and `fpm test` commands could not be executed. The source
layout follows FPM automatic target discovery, and `fpm.toml` was parsed as
valid TOML. All sources and targets were instead compiled directly with the
strict GNU Fortran command above.

## Scope of validation

The tests validate compilation, bounds and runtime checks, deterministic core
formulas, numerical sanity, output shapes, and bootstrap execution. They are
not a claim of exhaustive equivalence for every finite-sample edge case or of
bit-identical bootstrap output relative to R.
