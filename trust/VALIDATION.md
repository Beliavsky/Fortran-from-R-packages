# Validation

The translated source is validated with GNU Fortran using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. Rosenbrock minimization.
2. The rare hard-hard trust-region case from the upstream quartic saddle family.
3. Explicit easy-easy and hard-easy subproblem classification.
4. `parscale` behavior.
5. Maximization.
6. Restricted-domain/infeasible trial handling.
7. Iteration-limit behavior on an unbounded saddle.

The examples cover Rosenbrock minimization and concave-quadratic maximization.

FPM is not installed in the validation container, so the source/test/example tree is compiled directly with the strict flags above.  The `fpm.toml` manifest is parsed independently with a TOML parser, and the release archive is unpacked into a fresh directory and rebuilt from only its contents.
