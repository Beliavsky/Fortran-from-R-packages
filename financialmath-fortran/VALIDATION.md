# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64

FPM was not installed in the build environment. The manifest was parsed as
TOML and the project was compiled using the same automatic `src`, `app`,
`example`, and `test` layout expected by FPM.

## Compiler settings

Checked build:

```text
-std=f2018 -O0 -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace
```

Optimized build:

```text
-std=f2018 -O2 -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror
```

## Test suites

```text
test_cashflows: PASS
test_annuities: PASS
test_loans_bonds: PASS
test_derivatives: PASS
validation: PASS
```

The demo and both examples compile and execute under both builds.

## Reference checks

- Effective nominal-rate conversion is checked against direct compounding.
- TVM values use fixed analytical references.
- NPV and a 10 percent one-period IRR are checked directly.
- Cash-flow duration and modified-duration identities are checked.
- Level, arithmetic, and geometric annuities are checked against direct sums.
- Unknown payment, rate, increment, and growth parameters are recovered.
- Level, arithmetic, and geometric perpetuity identities are checked.
- A 30-year monthly amortization schedule pays the balance to zero and returns
  total principal equal to the original loan.
- A fixed coupon bond is checked against a closed-form reference price.
- Black-Scholes call and put values are checked against standard references.
- Put-call parity and call/put delta parity are checked.
- Long/short option tables and bull/bear spreads cancel exactly.
- Forward and prepaid-forward values are checked under no dividends and a
  continuous dividend yield.

## Release audits

- `fpm.toml` parses successfully as TOML.
- All translated text is ASCII.
- Maximum free-form Fortran line length is at most 132 columns.
- Every Fortran unit uses `implicit none` directly or through its program/module
  declaration scope.
- Every Fortran file has the preserved SPDX license and upstream attribution.
- Original and translated SHA-256 manifests are included.
