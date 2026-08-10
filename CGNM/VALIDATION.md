# Validation

The translated source is validated with GNU Fortran using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs cover:

1. algorithm-version-3 exact linear-model convergence;
2. legacy version-1 convergence;
3. finite-bound transformations and multi-objective targets;
4. low-level regularized/unregularized CGNR;
5. residual bootstrap and EBE extensions;
6. postprocessing helpers and R-type-7 column quantiles.

Both example programs are also compiled and run.

The source tree is scanned to ensure no translated free-form Fortran source
line exceeds 132 columns, and `fpm.toml` is parsed with Python's TOML parser.

The final release archive is extracted into a new directory and the same
strict build/test/example sequence is rerun solely from the archive contents.
