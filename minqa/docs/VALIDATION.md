# Validation

The test suite builds with Fortran 2018 conformance and runtime checking.

Covered cases:

- `newuoa`, `uobyqa`, and `bobyqa` recover the exact minimum of a shifted
  four-dimensional quadratic from the R package tests;
- BOBYQA minimizes the bounded two-dimensional Rosenbrock function;
- BOBYQA output remains inside every supplied bound;
- all three algorithms solve the six-variable Fletcher Chebyquad case used by the original package tests;
- the raw `390` maximum-evaluation result maps to public status `1`;
- invalid bounds are rejected before evaluating the objective;
- NaN objective values are converted to a large finite penalty;
- the example exercises all three public optimizers.

Build commands used by the supplied scripts:

```text
gfortran -std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace
gfortran -std=f2018 -O3
```

Warnings about exact floating-point comparisons are expected because those
comparisons are part of Powell's original trust-region algorithms and were
preserved intentionally.
