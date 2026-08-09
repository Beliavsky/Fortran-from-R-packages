# Translation coverage

## Translated

* R-facing QP model semantics: `P`, `c`, equality `A,b`, two-sided `G` bounds,
  and variable bounds.
* PIQP-style infeasible primal-dual predictor/corrector interior-point loop.
* Proximal primal and augmented-Lagrangian dual regularization (`rho`, `delta`)
  with adaptive reduction.
* Mehrotra affine predictor and centering/corrector step.
* Separate primal and dual fraction-to-boundary step lengths and `tau`.
* Dense regularized KKT elimination into an SPD solve.
* Cholesky factorization with iterative refinement.
* PIQP settings, result vectors, status values, objective/residual/gap fields.
* Persistent setup/solve/update model with fixed dimensions.
* Dense one-shot solve and Matrix-fortran CSC input adapter.
* Infinite constraint disabling and restoration of inactive slacks.
* Immediate detection of inconsistent two-sided/box bounds.

## Architectural differences from upstream PIQP 0.6.2

### Sparse backend

Upstream PIQP contains several specialized sparse KKT formulations and custom
sparse LDLT/orderings. This first Fortran release accepts CSC matrices through
the supplied Matrix-fortran package but converts them to dense arrays before
the KKT solve. Therefore sparse input is functionally supported, but very large
sparse problems do not yet obtain upstream PIQP's memory/time scaling.

### KKT implementation

The Fortran implementation eliminates equality variables and solves a
regularized SPD primal system by Cholesky. Upstream PIQP has dedicated KKT
classes, multiple sparse elimination strategies, factorization retries, and
more elaborate static regularization.

### Ruiz equilibration

The upstream dense/sparse preconditioners implement iterative Ruiz
preconditioning and optional cost scaling. The setting fields are retained for
API compatibility, but the first release does not yet reproduce the full Ruiz
scaling pipeline. The core regularization and iterative-refinement mechanisms
are active.

### Infeasibility certificates

Obvious inconsistent bound pairs are detected immediately. Upstream PIQP also
has proximal-progress based general primal/dual infeasibility detection. The
full certificate logic is not yet reproduced, so difficult infeasible problems
may end with `PIQP_MAX_ITER_REACHED` instead of `PIQP_PRIMAL_INFEASIBLE` or
`PIQP_DUAL_INFEASIBLE`.

### Warm starts on update

The model preserves data/settings across updates, but the current Fortran
`solve()` reinitializes primal-dual iterates rather than retaining all upstream
internal factorization/preconditioner/warm-start state. Numerical answers and
the public update semantics are preserved; allocation-free re-solves are not.

## Omitted R infrastructure

S7 classes, external pointers, `.Call`/RcppEigen registration, R coercion,
`simple_triplet_matrix` methods, formula-like R list handling, documentation
objects, and vignettes are not computational solver code and are not ported.
The original package is preserved under `original/`.
