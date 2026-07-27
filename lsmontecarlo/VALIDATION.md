# Validation record

The project was validated with GNU Fortran 14.2.0 in Fortran 2018 mode.

## Strict diagnostic build

The source, tests, application, and examples were compiled using:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace
-ffree-line-length-132
```

All tests passed:

```text
test_european: PASS
test_simulation: PASS
test_american: PASS
test_exotics_surface: PASS
```

## Test content

- Black-Scholes call and put values are checked against standard high-precision
  references.
- Put-call parity, zero-volatility pricing, and maturity payoffs are checked.
- Deterministic GBM paths are checked analytically.
- Antithetic path products are checked against their exact identity at every
  time step.
- Correlated GBM paths with `rho=1` are checked for exact equality.
- `first_value_row` is checked on a multirow matrix.
- Plain, antithetic, and control-variate American put prices are checked against
  a 1,200-step Cox-Ross-Rubinstein American put tree.
- The American price is checked against European and payoff bounds.
- A deterministic multiplier in the quanto model is checked against the
  corresponding scaled plain-American result. This is also a regression test
  for rank-deficient basis handling.
- Asian and quanto pricing, result metadata, and price surfaces receive bounds,
  shape, and summary checks.

The demo and both examples also compiled and ran successfully.

## Optimized build

The same source graph and executables were additionally compiled at `-O2` and
all tests passed.

## FPM

FPM was not installed in the validation environment. `fpm.toml` was parsed as
TOML and the project follows FPM automatic discovery for `src`, `app`,
`example`, and `test`. Direct GNU Fortran scripts are included for independent
validation.

## Release audits

The release tree was checked for:

- valid TOML syntax and the expected package/license fields
- ASCII-only translated text
- GPL SPDX identifiers and original attribution in every Fortran file
- `implicit none` in every Fortran compilation unit
- free-form Fortran lines no longer than 132 columns
- exact equality of the retained original tree to the supplied extraction
- successful verification of all provenance checksum manifests
- absence of object files, module files, executables, and build directories

The translated project contains 1,816 lines across library, application,
example, and test Fortran sources.

The final ZIP archive was extracted into a clean directory and the complete
strict validation script was rerun successfully from that extraction.
