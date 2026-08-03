# Porting notes

## Result and argument mapping

R lists are represented by `type(fingraph_result)`. R character initialization
choices are represented by the optional `initialization` argument, and an
explicit R numeric `w0` is represented by `initial_weights(:)`.

The R argument `d`, which may be scalar or vector, is represented by two
mutually exclusive optional arguments:

- `degree`: one common scalar degree
- `degrees(:)`: one degree per node

R dots were converted to underscores where necessary. The graph operator names
`L`, `A`, `D`, `Lstar`, `Astar`, and `Dstar` are valid Fortran identifiers and
were retained.

## Reused dependency code

Fingraph calls internal and exported routines from spectralGraphTopology and
uses fitHeavyTail in its documented workflows. To keep the FPM project
self-contained:

- Graph operators, graph metrics, dense pseudoinversion, symmetric Jacobi
  eigendecomposition, and nonnegative QP initialization were adapted from the
  completed `spectralGraphTopology-fortran` project.
- The deterministic RNG and multivariate Gaussian/Student-t simulation support
  used by tests and examples were adapted from `fitHeavyTail-fortran`.

The Fingraph estimators themselves remain direct translations of the three R
source files retained under `original/R/`.

## Eigenvalue ordering

R's `eigen(..., symmetric=TRUE)` returns eigenvalues in decreasing order. The
Fortran Jacobi routine returns them in increasing order. Index ranges were
translated accordingly:

- The `k` null-space vectors are columns `1:k` in Fortran.
- The positive-rank part of a `k`-component Laplacian is columns `k+1:p`.

This ordering change is essential to the k-component ADMM updates.

## Numerical safeguards

The translation preserves the algorithms while adding guards for conditions
that can otherwise create NaNs or invalid array sections:

- zero-row initial adjacency normalization
- nonpositive Student-t degrees of freedom
- constant columns during standardization
- logarithms of nonpositive eigenvalues
- exhaustion of `maxiter` without convergence
- zero Laplacian norm in relative-change tests

When `maxiter` is exhausted, the final iterate is returned with
`status=fg_no_convergence` and correctly sized histories.

## `record_objective`

The original R implementation accepts `record_objective` but computes and
stores the augmented Lagrangian unconditionally. The Fortran implementation
honors the documented argument: the k-component Lagrangian history is allocated
only when `record_objective=.true.`.

## Dense implementation

The original R packages use MASS and compiled linear algebra dependencies. This
translation intentionally uses dependency-free dense algorithms. It is suitable
for moderate node counts and verification work, not yet optimized for very
large sparse market graphs.
