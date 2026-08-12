# Validation

The translation was compiled with GNU Fortran using both a strict debug build
and an optimized build.

Strict flags:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

The sources fit within the standard free-form 132-column limit; no
`-ffree-line-length-none` option is required.

Regression reference values include:

- Steinhauser production LP: objective 93600, solution `(44,24,0)`, duals
  `(0,240,28.8)`, 0 Phase-I and 2 Phase-II iterations.
- Witte/Deppe/Born feed LP: objective 10.454545455, solution
  `(1.818181818,2.954545455)`, 2 Phase-I and 0 Phase-II iterations.
- Equality compatibility case: internal solver objective 1998 versus correct
  lpSolve-backed objective 1404 at `(41,33)`.
- Explicit dual solution of the production LP: `(0,240,28.8)` with both
  backends.
- MPS round-trips preserve objective, constraints, coefficients, and solution.
