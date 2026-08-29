# Validation

The automated test program covers:

- the original package's four-variable bounds example, recovering
  `[0, 0.75, 1.5, 2.25]`;
- the two-parameter Rosenbrock problem with analytic gradients;
- the same problem with bound-aware finite-difference gradients;
- Fletcher Chebyquad problems with 2, 3, 5, and 8 parameters, using the
  reference solutions in the package tests;
- the package's 100-parameter generalized Rosenbrock example;
- polymorphic callback data;
- evaluation limits, invalid bounds, non-finite callbacks, and monitor
  cancellation.

Both runtime-checked and optimized GNU Fortran builds are run by the release
validation scripts.

## Independent bounded Rosenbrock comparison

With two-sided bounds `[-2, 2]`, memory 5, `factr=1e7`, `pgtol=1e-9`, and
parameter-change stopping disabled, the Fortran example returns:

```text
x = 0.9999997588, 0.9999995521
f = 1.77654170e-13
iterations = 35
```

SciPy's L-BFGS-B implementation with the same settings returns the same
iteration count and agrees in the displayed parameters and objective value.
