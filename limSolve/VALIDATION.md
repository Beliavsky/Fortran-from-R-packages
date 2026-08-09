# Validation

The translated library, vendored dependencies, tests, and examples were built
with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

## Regression programs

Seven test programs pass:

1. `test_inverse` — generalized inverse, Lawson-Hanson-style NNLS, LDP.
2. `test_lsei` — unconstrained/equality LSEI, LDEI, covariance/rank output,
   and the upstream documented constrained LSEI example.
3. `test_linear` — tridiagonal and banded systems.
4. `test_block` — almost-block-diagonal reconstruction/solve.
5. `test_ranges` — LP, `xranges`, `varranges`, `varsample`.
6. `test_resolution` — row/column resolution and numerical rank.
7. `test_sampling` — mirror, random-direction, and coordinate-direction
   feasible-space sampling.

The upstream LSEI documentation example gives:

```text
x = 1.2424242424  3.0000000000  -0.7575757576
||A*x-B||^2 = 98.0606060606
```

## Examples

Both examples run successfully:

- `constrained_least_squares`
- `sample_polytope`

## Portability checks

- All translated and vendored library code compiles with
  `-Werror=implicit-interface`.
- No translated Fortran source line exceeds the standard free-form 132-column
  limit.
- The root and vendored `fpm.toml` files parse as TOML.
- FPM itself was not installed in the validation environment, so the identical
  FPM source/dependency tree was compiled directly with gfortran.
- The final archive was extracted to a new directory and the strict test script
  was rerun using only files contained in the archive.
