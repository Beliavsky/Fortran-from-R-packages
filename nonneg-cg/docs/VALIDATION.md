# Validation

## Direct C-kernel regression

The supplied `src/nonnegcg.c` kernel was compiled with a small standalone BLAS
shim and run on the package's Rosenbrock example from `x = (0, 2)`.

After exactly five iterations, the C kernel produced:

```text
status = 2
niter = 5
nfeval = 24
x = 1.1822794465797628, 1.4201196010889658
f = 0.083110622851909594
```

The Fortran translation reproduces these values to floating-point rounding.
The regression is deliberately stopped after five iterations because the
backtracking acceptance decision becomes sensitive to compiler-level rounding
near the final solution.

## Full convergence

The full Rosenbrock test verifies:

- tolerance status;
- feasibility throughout the solve;
- both variables within `1e-3` of one;
- final objective below `1e-8`;
- consistency between the legacy and actual objective-call counters.

## Additional tests

- quadratic optimum with one active lower bound;
- stationary point at the nonnegative boundary;
- optional polymorphic user data;
- maximum iteration and evaluation limits;
- monitor cancellation;
- invalid starting points;
- non-finite gradient detection.
