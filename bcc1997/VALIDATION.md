# Validation report

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- No external numerical libraries

## Checked build

The checked build used:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

## Optimized build

The complete suite was also rebuilt with `-O2`, `-Wall`, `-Wextra`,
`-Wimplicit-interface`, and `-Werror`.

## Test suites

```text
test_black_scholes_limit: PASS
test_pricing: PASS
test_strikes: PASS
test_validation: PASS
```

The demo and both examples compile and run in checked and optimized builds.

## Numerical references

Independent Python/SciPy integration of the original R formulas produced the
following references:

| Case | Call | Put |
|---|---:|---:|
| Upstream example 1 / Black-Scholes limit | 8.4333186900911 | 7.4383020650079 |
| Upstream example 2 | 8.6001374312277 | 7.6051208061445 |
| Fully active rate/variance/jump case | 7.5828223156054 | 9.3523119356215 |

The small difference between the first reference above and the analytical
Black-Scholes value is due to the upstream example's use of `1e-7` rather than
zero for two model volatilities. The Fortran result agrees with analytical
Black-Scholes pricing to better than `2e-11` in price.

Additional tests cover:

- characteristic-function values at a fixed transform argument;
- put-call parity;
- call/put monotonicity across strikes;
- jump-parameter irrelevance when jump intensity is zero;
- vector strike pricing;
- zero maturity and invalid parameters;
- fixed Black-Scholes values.

## FPM

FPM was not installed in the validation container. `fpm.toml` was parsed as
TOML, and the standard `src`, `app`, `example`, and `test` layout was compiled
directly with GNU Fortran. The project is intended to run with:

```text
fpm build
fpm test
fpm run bcc1997_demo
```
