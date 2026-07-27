# Validation record

This record is finalized by `scripts/build_validate.sh` and the release audit.

## Compiler configuration

- Compiler: GNU Fortran 14.2.0
- Language mode: Fortran 2018
- Flags: `-Wall -Wextra -Werror -fcheck=all -fbacktrace -O0`
- All modules, tests, the demo, and examples are linked from a clean build.
- Release Fortran line count: 3,227 across 17 source, test, demo, and example files.

## Test suites

- `test_distributions`: central/noncentral t and F reference values, inverse
  CDF identities, Sharpe and optimal-Sharpe distribution round trips.
- `test_estimation`: sample Sharpe ratios, higher moments, intervals, Markowitz
  weights, optimal Sharpe ratio, spanning delta, and SRIC.
- `test_hypothesis`: one-sample, paired/unpaired, equality, optimal, maximum,
  conditional, and power-test paths.
- `test_inference`: F/T-squared transformations, noncentrality estimators,
  optimal-Sharpe inference, and achieved-SNR confidence intervals.
- `test_unified`: empirical and Gaussian second-moment covariance estimators
  and inverse-moment delta propagation.

Expected successful output:

```text
test_distributions: PASS
test_estimation: PASS
test_hypothesis: PASS
test_inference: PASS
test_unified: PASS
```

## Independent distribution references

The deterministic constants in `test_distributions.f90` were independently
computed with SciPy and then tested in Fortran. Representative checks include:

```text
nct CDF(1.2; df=10, ncp=0.7) = 0.6751651489644870
nct PDF(1.2; df=10, ncp=0.7) = 0.3300988958641937
ncf CDF(2; df1=4, df2=20, ncp=3) = 0.6254295719853862
ncf PDF(2; df1=4, df2=20, ncp=3) = 0.2646774318185070
```

The test tolerances range from approximately `2e-12` to `3e-9`, depending on
whether the route uses direct evaluation or nested numerical inversion.

## Static release audits

The release audit verifies:

- `fpm.toml` parses as TOML.
- Every Fortran source contains `SPDX-License-Identifier:
  LGPL-3.0-or-later`.
- Every Fortran program unit uses `implicit none`.
- No translated Fortran line exceeds 132 columns.
- Release-authored text is ASCII.
- The original package and source archive checksums are retained.
- Supplied source archive SHA-256: `58ce955bd22c35caca3d42684aa9a67bbedc1af36f3b1875b6c3cc0db6386fd7`.
- No compiler objects, module files, executables, or transient build trees are
  included in the release archive.

## FPM note

The `fpm` executable was not installed in the translation environment. The
manifest was therefore parsed and audited, while compilation and execution
used the direct strict GNU Fortran build above. The project follows FPM's
automatic `src`, `app`, `example`, and `test` target layout.
