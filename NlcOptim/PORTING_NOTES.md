# Porting notes

## Numerical method

The upstream package follows Chapter 18 of Nocedal and Wright: it forms a
quadratic model of the Lagrangian, solves a constrained QP, performs a merit
line search, and updates the Hessian approximation. The Fortran translation
retains that architecture while making several numerical details explicit:

1. The QP is written in the supplied `quadprog` convention
   `A(:,j)^T d >= b(j)`.
2. Bounds and `<=` constraints are sign-transformed into that convention.
3. The BFGS update is Powell-damped so the QP matrix remains positive
   definite.
4. If the exact linearization is inconsistent, an elastic QP introduces one
   nonnegative maximum-violation variable. This replaces the upstream private
   `solqp` active-set fallback.
5. The merit function is an L1 exact penalty with Armijo backtracking.

## Compatibility choices

- Forward finite differences are the default, matching the upstream code.
- The starting point is projected into finite bounds before evaluation.
- Equality multipliers are sign-converted back to the conventional
  `f + lambda*h` Lagrangian.
- Constraint callbacks must return a fixed number of constraints throughout a
  solve. A changed size is reported as an error instead of relying on R's
  recycling behavior.
- Variables are represented as a one-dimensional vector. R's arbitrary input
  matrix shape was an interface convenience and is not part of the optimizer.

## Deliberate improvements

- Optional central differences and analytic derivative callbacks
- Explicit nonfinite-input/evaluation checks
- Structured status codes and messages
- Elastic feasibility restoration
- KKT and maximum-constraint diagnostics
- Reproducible checked and optimized build targets

## Dependency

The supplied `quadprog-fortran` translation is compiled from the three files
in `src/`. Its complete original translation is retained under `vendor/`, and
the uploaded archive is retained under `original/archives/`.
