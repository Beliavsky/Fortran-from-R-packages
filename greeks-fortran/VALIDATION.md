# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Python 3.13.5
- NumPy and SciPy 1.17.0 for independent analytical references

FPM itself was not installed in the validation container. The manifest was
parsed as TOML and the standard `src`, `test`, `app`, and `example` automatic
target layout was validated by direct compiler builds.

## Compiler configurations

Checked build:

```text
-std=f2018 -O0
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

Optimized build:

```text
-std=f2018 -O2
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

## Test results

```text
test_american: PASS
test_asian_implied: PASS
test_european: PASS
test_malliavin: PASS
```

The demo and both examples also compile and run in checked and optimized
configurations.

## Numerical checks

### European and binary options

- Fixed Black-Scholes put value, delta, gamma, and vega references
- Put-call parity with continuous dividend yield
- Asset-or-nothing minus strike times cash-or-nothing decomposition
- Binary delta checked by centered finite differences
- Vanilla vega checked by centered finite differences
- All higher-order analytical formulas compile and are available through the
  typed result

The fixed European put reference
`10.132463137873106` was independently reproduced with SciPy's normal CDF.

### Geometric Asian options

- Fixed put value, delta, and gamma references
- Vega checked by centered finite differences
- European and geometric-Asian implied-volatility recovery

The fixed geometric-Asian put reference
`6.566091953542269` was independently reproduced from the continuous
lognormal-average distribution with SciPy's normal CDF.

### American options

- A 500-step corrected American put value
- American put value no lower than the European put
- No-dividend American call equal to its European value
- Sensitivity sign checks
- American implied-volatility recovery

An independent NumPy CRR implementation produced:

```text
raw American tree:       6.088810110703037
raw European tree:       5.569527586515775
exact European value:    5.573526022256971
corrected American:      6.092808546444233
```

This matches the Fortran result.

### Malliavin Monte Carlo

- European call value and delta checked against exact Black-Scholes values
  using pathwise standard-error tolerances
- Geometric-Asian value checked against the exact analytical value
- Arithmetic-Asian value and direct/Malliavin delta smoke checks
- Regression control variate verified to reduce the fair-value standard error
- Deterministic seed behavior exercised
- Antithetic Brownian sampling exercised

Monte Carlo tests use statistical tolerances rather than machine-precision
fixed outputs so they remain portable across compilers and math libraries.

## Source and release audits

- 1,696 Fortran lines across 18 source, test, demo, and example files
- Every Fortran file contains `implicit none`
- Every Fortran file contains `SPDX-License-Identifier: MIT`
- All translated source and documentation are ASCII
- No translated free-form Fortran line exceeds 132 columns
- `fpm.toml` parses successfully
- Original source and translated-file SHA-256 manifests are included
- The final ZIP is extracted into a clean directory and the checked and
  optimized validation sequence is rerun before delivery
