# Translation coverage

## Computational behavior

| Upstream feature | Fortran 0.4.0 status |
|---|---|
| Dense ECOS model | Implemented |
| Native CSC/CSR storage | Implemented |
| Positive orthant | Implemented |
| Second-order cones | Implemented |
| Exponential cones | Implemented |
| Equality constraints | Implemented |
| Boolean/integer branching | Implemented |
| Setup/solve/update/cleanup | Implemented |
| Dense and sparse matrix updates | Implemented |
| Primal/dual/slack outputs | Implemented |
| Sparse symmetric KKT assembly | Implemented |
| Sparse symbolic/numeric LDL | Implemented |
| Predictor/corrector numeric-factor reuse | Implemented |
| Symbolic reuse over IPM iterations | Implemented |
| Symbolic reuse across workspace solves | Implemented |
| Same-pattern matrix-value symbolic reuse | Implemented |
| Warm-start reuse across workspace solves | Implemented |
| Iterative refinement | Implemented |
| Approximate minimum-degree ordering | Implemented |
| RCM ordering | Implemented as fallback |
| Cone-preserving equilibration | Implemented |
| Dynamic regularization | Implemented |
| Best-iterate restoration | Implemented |
| Inaccurate-optimal status | Implemented |
| Sparse primal-unbounded ray | Implemented |
| Sparse primal-infeasible dual ray | Implemented |
| Dual exponential-cone certificate transform | Implemented |
| Sparse ECOS_BB nodes | Implemented |
| ECOS_BB symbolic/warm reuse | Implemented |
| Sparse profiling counters | Implemented |
| Integrated ECOS homogeneous self-dual embedding | Not yet |
| Exact ECOS/Nesterov-Todd cone scaling | Not yet |
| Full SuiteSparse AMD implementation | Not yet |

## What changed from v0.2

### Ordering

v0.2 used reverse Cuthill-McKee. v0.4 defaults to a native approximate-minimum-degree ordering with sparse dynamic adjacency lists and lazy heap updates. Exact clique fill is tracked for ordinary fronts; large fronts use bounded-cost approximate degree updates. RCM remains available.

The implementation is intentionally described as **AMD-style**, not as a direct SuiteSparse AMD port. SuiteSparse AMD has additional supervariable, quotient-graph, aggressive-absorption, and dense-row machinery that is not reproduced here.

### Equilibration

Sparse problems now receive iterative scaling before the IPM solve. Variable columns are scaled individually. Equality rows and linear-cone rows can scale independently. Each SOC block and each exponential-cone triple is scaled by one positive scalar, preserving cone geometry.

Solutions and duals are mapped back to original coordinates. This is a Ruiz-style conditioning pass, not a claim of bit-identical ECOS preprocessing constants.

### Persistent cache

`ecos_workspace` owns a sparse cache containing:

- permutation;
- elimination tree;
- symbolic `L` column pointers/capacity;
- KKT structure hash;
- previous scaled `x`, `y`, scalar slacks, and scalar duals.

`c/h/b` changes preserve the cache. Same-pattern matrix value changes preserve symbolic data but invalidate the warm point. Structural changes invalidate symbolic data.

### Sparse certificates

Large sparse failure diagnostics no longer require dense conversion. Two homogeneous feasibility problems are used:

1. primal ray for unboundedness: `A*d=0`, `-G*d in K`, `c^T*d=-1`;
2. dual ray for primal infeasibility: `A^T*y+G^T*z=0`, `b^T*y+h^T*z=-1`, `z in K*`.

For exponential dual blocks `(u,v,w)`, the certificate model uses `(-v,-u,e*w) in K_exp`.

These are normalized sparse certificates, but they are solved **after** an unsuccessful main solve. Upstream ECOS obtains certificates naturally from its integrated homogeneous self-dual embedding.

### Sparse ECOS_BB

Sparse integer problems remain sparse at every node. Two structurally fixed bound rows are reserved per integer variable. Inactive rows contain explicit structural zeros, so branch decisions change values rather than sparsity. This enables symbolic factorization reuse across nodes.

## Remaining differences from upstream ECOS

### Homogeneous self-dual embedding

The main Fortran IPM remains a conventional primal-dual predictor/corrector method. v0.4 adds sparse homogeneous ray-certificate problems, but it does not augment every Newton system with ECOS's `tau/kappa` homogeneous self-dual variables.

Consequences:

- optimal trajectories/iteration counts need not match ECOS;
- infeasibility certificates are post-failure auxiliary solves rather than native HSD iterates;
- not every upstream inaccurate primal/dual infeasibility status is reproduced.

### Cone scaling

SOC and exponential cones have specialized sparse derivatives/curvature, but the solver still represents them through smooth scalar convex inequalities. It does not yet reproduce ECOS's exact Nesterov-Todd SOC scaling blocks or its precise exponential-cone scaling/KKT layout.

### AMD

The default ordering is a native approximate-minimum-degree implementation, not SuiteSparse AMD itself. It is substantially better than identity/RCM on the irregular fill regression used here, but fill counts can differ from ECOS.

### Sparse KKT numeric updates

The elimination tree/pattern is reused, but numeric KKT values are rebuilt into sparse triplet/CSC storage each Newton iteration. Upstream ECOS has more specialized preallocated update paths. Eliminating those remaining allocations is future performance work.

### Regularization/refinement

v0.4 uses signed quasi-definite pivot regularization with a dynamic retry schedule and residual-improvement iterative refinement. Constants and all fallback branches are not bit-identical to ECOS.

### ECOS_BB

Sparse node reuse is implemented, but original ECOS_BB node heuristics, rounding heuristics, pseudocost-style choices, and every incumbent/cutoff detail are not reproduced exactly.

### Small dense diagnostic fallback

After sparse certificate attempts, problems smaller than `dense_diagnostic_limit` may still use the legacy dense diagnostic fallback. Set `dense_diagnostic_limit=0` in `ecos_settings` to disable it completely.

## MatrixExtra

The optional MatrixExtra adapter is included for sparse construction/interoperability. MatrixExtra's current translated direct sparse-solvers convert to dense internally, so they are deliberately not used for ECOS KKT factorization.

## R-only functionality omitted

- `.Call`/SEXP interfaces;
- R external pointers/finalizers;
- S3 printing;
- R `Matrix`/`slam` dispatch;
- R data-frame status descriptions.
