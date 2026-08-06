# Build report

Commands used for validation:

```text
make check
make optimized
make example
```

The checked configuration uses Fortran 2018, `-Wall -Wextra -Werror`, runtime
checking, backtraces, and explicit-interface enforcement. The optimized
configuration uses `-O3` with the same warning policy.

FPM was not installed in the validation environment. `fpm.toml` was parsed as
TOML and the same source dependency order was compiled through the Makefile.
