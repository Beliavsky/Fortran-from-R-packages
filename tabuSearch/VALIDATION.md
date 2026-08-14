# Validation

The release was compiled with GNU Fortran 14.2.0 using

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

All five permanent test programs pass:

1. `test_reference` checks a hand-derived deterministic preliminary-search
   trace and exact optimum for a six-bit weighted objective.
2. `test_negative_neigh` uses an all-nonpositive objective with `neigh < size`
   and verifies that every stored utility is exactly the objective of the
   stored binary configuration.
3. `test_summary` independently recomputes summary counts and move
   frequencies.
4. `test_repeat_rng` checks Park-Miller reference outputs, `repeat_all`
   history size, and seeded reproducibility.
5. `test_linear_stress` solves 100 independently generated 20-bit linear
   binary objectives.  Their exact global solution is known coordinatewise;
   all 100 searches found the exact optimum value and configuration.

The example also runs under the same checked build and returns the expected
best configuration for its weighted binary objective.

The FPM executable was not installed in the validation environment.  The
exact FPM source/test/example tree was therefore compiled directly with
`gfortran`; `fpm.toml` is parsed independently during the release audit.
