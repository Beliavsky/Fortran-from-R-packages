# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran standard mode: Fortran 2018
- Source layout: standard FPM automatic discovery
- Fortran source, tests, application, and examples: 1736 lines

The FPM executable was not installed in the validation environment. The
manifest was parsed independently as TOML, and all targets were compiled using
the same source separation used by FPM.

## Checked build

The checked build used:

```text
-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace -O0 -g
```

- warnings promoted to errors;
- bounds, allocation, pointer, and array-temporary runtime checking;
- backtraces;
- interface and conversion diagnostics; and
- no optimization.

Results:

```text
test_benders: PASS
test_projection: PASS
test_risk: PASS
test_simplex: PASS
validation: PASS
```

The demo and both examples also compiled and ran successfully.

## Optimized build

The complete source, tests, demo, and examples were rebuilt with `-O2` while
retaining strict warning and interface diagnostics. All four tests passed.

## Independent numerical references

### Risk statistics

For losses `(2, 3, 1, 1)` with equal probability and `alpha = 0.90`, the tests
verify:

- VaR = 3;
- CVaR = 3;
- mean = 1.75; and
- MAD = 0.75.

### Linear programming

The simplex suite verifies a bounded LP with optimum `(2, 2)` and objective
`-10`. It also verifies detection of an infeasible problem containing both
`x <= 0` and `x >= 1`, exercising Phase-I artificial variables.

### Benders decomposition

A five-scenario, two-asset problem was evaluated independently as a full linear
program. The tests verify:

| Risk | First weight | Optimized risk |
|---|---:|---:|
| CVaR | 0.2631578947368421 | 0.0115789473684211 |
| deviation CVaR | 0.2631578947368421 | 0.0286315789473684 |
| LSAD | 0.0322580645161290 | 0.00993548387096774 |
| MAD | 0.0322580645161290 | 0.0198709677419355 |

The values agree with independent SciPy/HiGHS full-LP calculations. SciPy is
used only for offline reference generation and is not a project dependency.

### Optimal-face projection

The projection suite verifies:

- a unique MAD-optimal portfolio against an independent full-LP reference;
- budget, bound, and expected-return constraints;
- a degenerate zero-risk optimal face, where the benchmark portfolio is
  recovered; and
- a direct Zhao-Li compatibility-routine feasibility and residual check.

## Source audits

The release is checked for:

- valid TOML;
- ASCII-only translated text;
- free-form source lines no longer than 132 columns;
- `implicit none` in every Fortran program and module;
- GPL-3.0-only SPDX headers in every translated Fortran file;
- retained unmodified upstream source and source archive; and
- SHA-256 manifests for original and translated files.

## Reproduction

With FPM:

```text
fpm build
fpm test
fpm run portfoliooptim_demo
fpm run --example benders_risk_measures
fpm run --example benchmark_projection
```

Without FPM:

```text
sh scripts/validate.sh
```

On Windows:

```bat
scripts\validate.bat
```

The final ZIP archive was extracted into a clean directory and the checked and
optimized validations were repeated successfully.
