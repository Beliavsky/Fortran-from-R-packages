# Validation record

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Direct compilation because the FPM executable was not installed in the validation runtime
- FPM manifest parsed independently as TOML

## Checked build

Compiler options:

```text
-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace
```

Results:

```text
test_fixed_income: PASS
test_market: PASS
test_options: PASS
test_portfolio: PASS
test_processes: PASS
```

The demo and all three examples compile and execute.

## Optimized build

Compiler options:

```text
-std=f2018 -O2 -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fbacktrace
```

All five tests, the demo, and all three examples pass.

## Test coverage

- Generalized Black-Scholes fixed references, put-call parity, CRR convergence, tree consistency, Kirk spread references, and barrier expiry/knockout behavior.
- Deterministic GBM and OU paths, time-varying OU interpolation, exact discrete OU fitting, multivariate moments and covariance, and zero-intensity OU-jump behavior.
- Par bond pricing and duration, NPV terminal conventions, semiannual IRS schedule generation, commodity business-day weights, and direct/calendar swap prices.
- Simplex fixed optimum, refinery objective consistency, random portfolio weight constraints, and minimum-risk/maximum-Sharpe indexing.
- Absolute and relative returns, roll-day masks, full/bull/bear betas, performance statistics, and moving-average strategy accounting.

## Release size

The release contains 2,703 lines across 19 Fortran source, test, demo, and example files.

## Source audits

- All translated source and documentation text is ASCII.
- Every Fortran source has `implicit none` and an MIT SPDX header.
- Free-form Fortran lines are at most 132 columns.
- `fpm.toml` parses successfully and declares the MIT license.
- GARCH does not appear in translated Fortran source; the upstream `garch.R` remains only in the retained original package.
- Original and translated file checksum manifests are included under `provenance`.
