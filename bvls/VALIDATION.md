# Validation

The translated source is intended to compile with strict GNU Fortran checks:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Validation includes analytic bounded least-squares cases, infinite bounds,
warm starts, invalid/fixed bounds, and fixed numerical references obtained from
the original Stark-Parker fixed-form implementation.

During translation development, the modern solver was linked alongside the
original fixed-form source and compared on 1000 deterministic problems, both
cold and warm started. All coefficient vectors and residual norms agreed within
`5e-11`; the tested cases produced zero mismatches.
