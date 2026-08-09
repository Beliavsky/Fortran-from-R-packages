# Validation

Validation compiler: GNU Fortran 14.2.0.

Flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. primal cone/isotonic projection (`cone_a`),
2. dual/constraint-cone projection (`cone_b`),
3. ordinary and shifted-constraint quadratic programming,
4. shape-edge generation and monotonic projection,
5. constrained regression,
6. shape-restricted regression,
7. beta and Student-t numerical utilities,
8. QR/rank and irreducibility utilities.

Both examples are also compiled and run. The final release archive is re-extracted into a fresh directory and the same strict suite is rerun from only those contents.
