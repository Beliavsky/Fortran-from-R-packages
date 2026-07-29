# Validation report

## Toolchain

The project was compiled with GNU Fortran 14.2.0 in two configurations.

Checked build:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

Optimized build:

```text
-std=f2018 -O2
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

## Test suites

```text
test_math: PASS
test_portfolio: PASS
test_naive: PASS
test_sis: PASS
test_objectives: PASS
test_gh: PASS
```

The demo and both examples also compile and run under checked and optimized
builds.

## Coverage of tests

- Student-t CDF/quantile inversion.
- Gamma CDF/quantile inversion.
- Orthonormal-basis construction.
- Minimax allocation symmetry.
- Exact center-state t-copula portfolio return.
- Tail and excess response definitions.
- Rare-event boundary and direction optimization.
- Seeded naive-simulation reproducibility.
- Threshold monotonicity.
- Positive estimated variances and conditional excess.
- SIS versus high-sample naive Monte Carlo within combined simulation error.
- Multi-threshold conditional-excess optimization.
- Maximum-relative-error allocation.
- GH inverse-table monotonicity and seeded GH simulation.
- Runtime bounds and allocation checks.

## Archive validation

The release validation script rebuilds the source from the exact ZIP, reruns all
six tests, the demo, and both examples, and verifies the original and translated
SHA-256 manifests.

FPM was not installed in the construction environment. The manifest was parsed
as TOML and the project was validated using the same automatic `src`, `test`,
`app`, and `example` directory layout used by FPM.
