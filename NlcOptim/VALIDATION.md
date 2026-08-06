# Validation

Five deterministic test programs cover:

1. Unconstrained quadratic optimization with an analytic gradient.
2. Linear equality constraints combined with lower and upper bounds.
3. A nonlinear circle equality with an analytic Jacobian.
4. Nonlinear inequalities and Rosenbrock optimization from an infeasible start.
5. Direct validation of the vendored `quadprog` solver and invalid-input
   handling.

The demonstration translates the first nonlinear equality-constrained example
from the upstream documentation. In the validated GNU Fortran build it reaches
an objective near `5.395e-2` with maximum equality violation below `1e-6`.

Validated configurations:

- GNU Fortran 14.2, `-O0`, warnings as errors, runtime checks, backtraces, and
  floating-point traps.
- GNU Fortran 14.2, `-O3`, warnings as errors.

The exact iteration path need not match R because line-search safeguards,
finite-difference ordering, and the QP fallback differ. Solutions and KKT
conditions are the validation targets.
