# Validation

The translated library and tests were built with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Tests cover:

1. default All To One optimization of a convex quadratic;
2. T3A optimization;
3. Pareto optimization;
4. zero-migration initialization and result bookkeeping;
5. invalid option handling.

The upstream R package contains one tinytest checking only the default SOMA
strategy on `sum(x^2)`. The Fortran test suite therefore extends coverage to
all exported computational strategies.

Both included examples are also compiled and run under the same flags.

FPM was not installed in the validation container, so `fpm build` itself could
not be executed there. The FPM source layout and `fpm.toml` are supplied, and
the same source tree is compiled directly by the strict validation scripts.
