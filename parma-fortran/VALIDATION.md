# Validation

## Environment

The release was validated on July 25, 2026 with:

```text
GNU Fortran (Debian 14.2.0-19) 14.2.0
```

FPM was not installed in the validation environment. The `fpm.toml` file was
parsed with Python's TOML parser and follows FPM automatic discovery for the
`src`, `app`, `example`, and `test` directories. Compilation was performed
manually using the same source graph declared by that layout.

## Compiler checks

Every library module, test, application, and example was compiled with:

```text
-std=f2018
-pedantic
-Wall
-Wextra
-Werror
-Wimplicit-interface
-Wconversion-extra
-Wcharacter-truncation
-fcheck=all
-fbacktrace
```

The project contains 24 Fortran files and 3,582 lines of Fortran source across
`src`, `app`, `example`, and `test`.

GNU `ld` reports that objects using nested Fortran procedure callbacks require
an executable stack. This is caused by GNU Fortran's trampoline implementation
for internal procedures. It did not prevent linking or execution in the test
environment, but hardened systems that prohibit executable stacks may need the
callbacks moved to non-nested module procedures.

## Test results

```text
test_core: PASS
test_optimizers: PASS
test_socp: PASS
test_features: PASS
```

The tests cover deterministic risk formulas, entropy, quantiles, lower partial
moments, lag/date helpers, CARA gradients, QP, LP, binary MILP, CMA-ES, SOCP,
high-level portfolio optimization, threshold-999 behavior, string-based risk
dispatch, benchmark covariance, covariance-only optimization, leverage,
cardinality, and simulated feasible weights.

The following executable targets also compiled and ran successfully:

```text
parma_demo
efficient_frontier
risk_measures
```

## Release audits

The following release checks passed:

```text
fpm.toml TOML parse: PASS
ASCII-only text files: PASS
GPL-3.0-or-later SPDX headers: PASS
Fortran line length <= 132: PASS
implicit none in every Fortran source: PASS
original archive SHA-256 match: PASS
original-tree checksum manifest: PASS
translated-file checksum manifests: PASS
clean release tree with no objects or build directories: PASS
archive extraction and strict rebuild: PASS
```

The exact release archive was extracted into a clean directory and all four
tests, the application, and both examples were rebuilt and executed again.

## Interpretation

The validation establishes that the translated code is internally consistent,
compiles under strict GNU Fortran 2018 diagnostics, and passes the included
reference and smoke tests. It does not claim bit-for-bit agreement with every
external R solver used by the original package. Nonunique optima and stochastic
algorithms can differ while satisfying the same mathematical formulation.
