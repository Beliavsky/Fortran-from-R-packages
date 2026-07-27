# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64
- FPM was not installed in the validation environment, so the included direct
  build scripts were used. `fpm.toml` was independently parsed as TOML and uses
  FPM automatic source, executable, example, and test discovery.

## Checked build

The checked build uses:

```text
-std=f2018 -O0 -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace
```

Results:

```text
test_covariance: PASS
test_formulas: PASS
test_inference_frontier: PASS
test_means: PASS
test_portfolios: PASS
test_random: PASS
validation (debug): PASS
```

## Optimized build

The complete source, tests, demo, and examples were rebuilt with `-O2` and the
same compile-time diagnostics:

```text
test_covariance: PASS
test_formulas: PASS
test_inference_frontier: PASS
test_means: PASS
test_portfolios: PASS
test_random: PASS
validation (optimized): PASS
```

## Numerical coverage

The tests include:

- Fixed sample covariance, BGP14 covariance shrinkage, and BGP16 inverse
  covariance shrinkage references.
- Positive-definiteness and symmetry of nonlinear Ledoit-Wolf shrinkage.
- Fixed Bayes-Stein, James-Stein, and BOP19 mean-shrinkage references.
- Fixed traditional MV, BDOPS21 MV-shrinkage, and BDPS19 GMV-shrinkage results.
- Weight-sum identities and `p > n` Moore-Penrose portfolio paths.
- Asymptotic MVSP testing and Bayesian frontier monotonicity.
- Seed reproducibility and prescribed eigenvalues for random covariance
  generation.
- Deterministic alpha and alpha-variance formulas.

Fixed reference values were independently calculated with NumPy using the
published formulas.

## Source audit

- 1,570 lines of Fortran across source, tests, demo, and examples.
- Every translated Fortran unit includes `implicit none`.
- Every translated Fortran unit includes `SPDX-License-Identifier:
  GPL-3.0-only`.
- Translated source is ASCII-only.
- Free-form lines do not exceed 132 columns.
- `fpm.toml` parses successfully.
- Original and translated file checksums are recorded under `provenance`.
