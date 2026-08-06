# Validation

The automated test program checks:

1. A 20-dimensional quadratic with all four exposed line-search modes.
2. The Rosenbrock function with the default More-Thuente search.
3. OWL-QN against an exact soft-threshold solution.
4. Partial OWL-QN regularization with unpenalized variables.
5. Separate objective and gradient callbacks.
6. Polymorphic user-data callbacks.
7. Progress-callback cancellation.
8. Invalid line-search and OWL-QN index handling.
9. Detection of an already minimized initial point.

With the package defaults, the Rosenbrock example starting at `(-1.2, 1)`
produces:

```text
x = 1.0000005950578863, 1.0000011922387717
f = 3.5454445020444373e-13
iterations = 37
evaluations = 45
```

This agrees with the numerical result printed in the supplied R package
vignette (`3.545445e-13`, parameters approximately `1.000001`).

The release is built in two configurations:

- runtime checks, bounds checks, and floating-point traps
- optimized `-O3` build
