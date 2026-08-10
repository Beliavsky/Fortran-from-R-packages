# Translation notes

## Scope

This port targets the **computational code** of R package `scs` 3.2.7 and the
SCS 3.x source bundled inside it. R object dispatch, `.Call` registration,
Matrix/slam coercion code, documentation generation, and test harness glue are
not reproduced in Fortran.

The public API is intentionally Fortran-native rather than an emulation of R's
lists and S3/S4 classes.

## Representation changes

### CSC sparse matrices

The C/R implementation ultimately presents SCS with compressed sparse column
arrays. `type(scs_matrix)` preserves that format but uses 1-based indices:

- `p(j)` is the first stored entry of column `j`;
- `p(n+1) = nnz + 1`;
- `i(k)` is the 1-based row of stored value `x(k)`.

### Positive-semidefinite cones

SCS's packed lower-triangular `svec` convention is preserved. Diagonal entries
are stored directly; off-diagonal entries are scaled by `sqrt(2)`. Projection
unpacks the symmetric matrix, performs a symmetric eigendecomposition, clips
negative eigenvalues, and repacks it.

### Box cone

The first entry is the box scaling variable `t`; the remaining entries satisfy
`t*bl <= x <= t*bu`, matching SCS's box-cone formulation.

## Solver translation

`scs_solver` follows the upstream homogeneous self-dual Douglas-Rachford flow:

1. validate input and cone dimensions;
2. optionally normalize `A`, `P`, `b`, and `c`;
3. form and factor the quasi-definite linear system;
4. perform the affine projection;
5. project onto the dual cone;
6. update the Douglas-Rachford iterate with over-relaxation;
7. periodically compute residuals and certificates;
8. optionally adapt scaling and use Anderson acceleration;
9. classify the final solution/certificate and unscale it.

Residuals are evaluated directly in the original problem coordinates. This
avoids exposing the C implementation's internal normalized-residual bookkeeping
while preserving the stopping criteria in original units.

## Linear-system backend in v0.2.0

v0.1.0 used a dense native LDL^T factorization of the SCS KKT matrix. v0.2.0
replaces that implementation with a sparse backend split into two modules:

- `scs_ldlt`: sparse KKT assembly, diagonal-update bookkeeping, symbolic-state
  ownership, and the SCS-facing factor/solve object;
- `scs_qdldl`: native Fortran translation of QDLDL's elimination tree, sparse
  numeric LDL^T factorization, and triangular solves.

The KKT matrix is the upper triangular part of

```text
[ P + R_x    A^T ]
[ A         -R_y ]
```

and is assembled directly into 1-based CSC storage. The diagonal is always
explicit, as required by QDLDL.

### Symbolic/numeric split

On the first `factorize` call, `scs_ldlt`:

1. assembles the sparse KKT pattern and values;
2. computes the QDLDL elimination tree and column counts;
3. allocates sparse `L` and numeric work arrays;
4. performs numeric factorization.

During SCS adaptive scaling, only `R_x`/`R_y` changes. Subsequent factorization
calls therefore update the stored KKT diagonal and rerun only the numeric QDLDL
phase. `test_sparse_ldlt` verifies that the symbolic-analysis counter remains
one while the numeric-factorization counter increments.

### Ordering

Upstream SCS first applies SuiteSparse AMD to the KKT sparsity graph. v0.2.0
uses natural KKT ordering. This was chosen deliberately to make the sparse
factorization available without prematurely translating the roughly 4,000-line
bundled AMD implementation.

Consequences relative to v0.1.0:

- dense `(n+m)^2` KKT/factor storage is eliminated;
- arithmetic is sparsity/fill dependent instead of fixed cubic dense work;
- adaptive rescaling reuses symbolic factorization state;
- large sparse problems with favorable natural ordering improve dramatically;
- problems with poor natural ordering may still create much more fill than
  upstream AMD+QDLDL.

Consequences relative to upstream SCS:

- QDLDL numerical semantics and data structures are closely matched;
- AMD fill reduction is still missing, so this release does **not** claim full
  upstream sparse-direct performance parity.

## QDLDL translation details

`scs_qdldl.f90` is translated from the QDLDL source bundled with SCS and keeps
its Apache-2.0 provenance. The main adaptation is systematic conversion from
0-based C indices to 1-based Fortran CSC indices. QDLDL's `-1` unknown-tree
sentinel becomes `0` because zero is outside valid Fortran matrix row indices.

The numeric routine additionally accumulates duplicate input entries into the
working row value, which is safe for canonical SCS input and more robust to
uncollapsed duplicates.

## Cone projections

`scs_cones` translates:

- zero and nonnegative projections;
- box-cone scalar Newton iteration;
- SOC projection;
- PSD projection;
- primal and dual exponential-cone projection;
- primal and dual power-cone projection;
- Moreau construction of the dual-cone projection.

The PSD eigensolver is a portable Jacobi method rather than BLAS/LAPACK so the
FPM package has no required dependency. This is a performance choice, not a
change to the cone definition.

## Anderson acceleration

`scs_acceleration` implements the limited-memory Anderson scheme used by SCS,
including positive-lookback/type-I and negative-lookback/type-II modes,
regularization, safeguard rejection, and history reset after a rescaling. The
small least-squares systems are solved internally with partial pivoting.

## Deliberately omitted facilities

- R `.Call` registration and R list construction;
- conversion from `Matrix`, `slam`, and ordinary R matrices;
- SuiteSparse AMD ordering (planned sparse-backend follow-up work);
- optional MKL and GPU backends;
- optional iterative/indirect linear solver;
- Ctrl-C signal-handler integration;
- CSV logging and problem serialization;
- build-system code for CMake/GNU Make;
- plotting (none of computational relevance in the package).

## Validation performed for v0.2.0

The complete solver regression suite is run with GNU Fortran 14.2.0 at `-O0`
with `-fcheck=all -fbacktrace`, and again at `-O2`.

The semidefinite regression uses the exact matrix, vectors, cone sizes, and
reference solution from `inst/tinytest/test_psd.R`.

The direct sparse backend test independently checks KKT solve residuals and
symbolic reuse. An included benchmark compares sparse QDLDL with the v0.1.0
dense LDL^T algorithm on a synthetic banded KKT system.
