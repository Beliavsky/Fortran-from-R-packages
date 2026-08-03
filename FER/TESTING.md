# Testing

The test target checks:

- Bachelier and Black-Scholes put-call parity.
- Recovery of known Bachelier and Black-Scholes volatilities.
- CEV limiting and probability behavior.
- SABR and NSVh reference cases.
- Margrabe and spread-option formulas.
- Finite and nonnegative outputs for representative parameter sets.

Run with:

```text
fpm test
```

For strict direct compilation with GNU Fortran, recommended flags are:

```text
-std=f2018 -pedantic -Wall -Wextra -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
-fimplicit-none -O2
```
