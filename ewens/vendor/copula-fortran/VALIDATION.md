# Validation

## Toolchain

The release was compiled with GNU Fortran 14.2.0 in Fortran 2018 mode.

Checked flags:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

Optimized flags:

```text
-std=f2018 -O2
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

## Test programs

```text
test_compositions_special: PASS
test_elliptical: PASS
test_families: PASS
test_fitting: PASS
test_simulation_empirical: PASS
```

The demo and both examples also compile and run under the checked build.

## Numerical references

Validation includes:

- Independent bivariate-normal CDF and copula-density reference values.
- Independent bivariate-Student CDF and copula-density reference values.
- Closed-form Clayton, Gumbel, Frank, FGM, Plackett, Marshall-Olkin, and
  rotation identities.
- Closed-form Kendall tau, Spearman rho, and tail-dependence identities.
- Three-dimensional Gaussian probability and seeded simulation smoke tests.
- Seeded Clayton sample Kendall tau.
- Rosenblatt/inverse-Rosenblatt round trips.
- Pseudo-observation tie handling and empirical-copula bounds.
- Deterministic permutation-test smoke checks.
- Inversion-of-tau and maximum-pseudo-likelihood fitting on seeded samples.
- Mixture, Khoudraji, and nested-Clayton identities.
- Stirling, Eulerian, Sibuya, and logarithmic-series identities.

## Source audit

The translated source is ASCII-only, uses free-form lines no longer than 132
columns, includes `implicit none`, and carries GPL-3.0-or-later SPDX headers.
The FPM manifest parses as TOML.

## FPM

FPM was not installed in the validation container. The manifest uses standard
automatic discovery for `src`, `test`, `app`, and `example`. Direct validation
scripts reproduce the complete build and test sequence on Linux and Windows.
