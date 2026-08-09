# Validation

The release is validated with GNU Fortran using:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs cover:

1. the upstream `test-solve-simple-qp.R` dense solution and multipliers;
2. the upstream update/re-solve values;
3. equality-only and two-sided inequality QPs;
4. a zero-Hessian linear objective with box constraints;
5. Matrix-fortran CSC input;
6. invalid settings and inconsistent bounds.

Both examples are compiled and run under the same flags. The strict script also
compiles every Matrix-fortran dependency module needed by the FPM target.

The final release archive is additionally extracted to a fresh directory and
rebuilt from only the archive contents before publication.
