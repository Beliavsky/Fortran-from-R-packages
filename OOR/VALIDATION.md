# Validation

The translated source was compiled with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs cover:

1. all exported benchmark functions;
2. deterministic 1-D SOO on the upstream guirland example;
3. deterministic 2-D SOO and retained tree/history output;
4. stochastic StoSOO with repeated leaf sampling;
5. POO with a reproducible scalar maximization problem.

Both examples are also compiled and executed with the same strict flags.

FPM was not available in the validation container, so the FPM project layout
and manifest are included but the numerical validation is performed by the
compiler-only scripts in `scripts/`.
