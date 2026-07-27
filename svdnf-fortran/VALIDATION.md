# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64

FPM was not installed in the translation environment. The FPM manifest was
validated as TOML and the source tree follows standard `src`, `app`, `example`,
and `test` discovery conventions. Direct GNU Fortran scripts are included.

## Checked build flags

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

## Optimized build flags

```text
-O2
-std=f2018
-Wall
-Wextra
-Wconversion-extra
-Wimplicit-interface
-Werror
```

## Test suites

```text
test_probabilities: PASS
test_models: PASS
test_filter: PASS
test_jumps: PASS
test_forecast: PASS
test_optimization: PASS
```

The demo and all three examples compile and run in checked and optimized builds.

## Tested behavior

- fixed normal, gamma, Poisson, and binomial references;
- Heston drift and diffusion equations;
- parameter packing and unpacking;
- every built-in grid family;
- custom callback simulation, filtering, and estimation;
- filter normalization at every observation;
- positive finite likelihood contributions;
- filtering and prediction percentiles;
- CAPM-SV factor adjustment;
- Bates and Pitt-Malik-Doucet jumps;
- Duffie-Pan-Singleton jump-size grids;
- normalized transition matrices;
- deterministic seeded forecasts;
- confidence interval ordering;
- built-in maximum-likelihood improvement;
- numerical Hessian dimensions and exact symmetry.

## Clean release validation

The final ZIP is extracted into a fresh directory and `scripts/validate.sh` is
run there. Source audits verify:

- ASCII-only translated text;
- free-form line lengths no greater than 132 columns;
- `implicit none` in every Fortran program or module;
- GPL-3.0-only SPDX headers;
- matching original and translated SHA-256 manifests.
