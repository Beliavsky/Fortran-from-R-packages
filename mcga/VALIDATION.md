# Validation

The translation was compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Tests cover:

1. bit-preserving scalar/vector double-byte round trips;
2. one- and two-point byte crossover and byte mutation wraparound;
3. the full mcga-specific real/byte crossover operator library;
4. single-objective MCGA minimization;
5. the upstream multi-objective rank-score formula;
6. multi-objective MCGA execution and final rank ordering.

Both examples are also compiled and executed by the strict validation script.

Some GNU/Linux linkers print a non-fatal executable-stack note when an example
or test passes an internal Fortran procedure as a callback. This is generated
by gfortran's trampoline implementation for contained procedure callbacks; the
library source itself compiles without warnings under `-Werror` and all
callback dummy procedures have explicit abstract interfaces.

FPM is not installed in the validation container. `fpm.toml` is parsed with a
TOML parser and the exact source/test/example tree is compiled independently.
The release archive is additionally extracted into a clean directory and the
strict script is rerun there.
