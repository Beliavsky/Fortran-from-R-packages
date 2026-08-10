# ECOSolveR-fortran

Modern Fortran/FPM translation of the computational surface of **ECOSolveR 0.6.1** and the ECOS conic model it wraps.

Version **0.4.0** builds on the sparse v0.2 backend with sparse ordering, cone-preserving equilibration, persistent symbolic/warm-start reuse, sparse infeasibility/unboundedness certificates, best-iterate/inaccurate-status handling, dynamic KKT regularization, profiling counters, and a sparse ECOS_BB path.

## Implemented

- `min c^T x` subject to `A x = b` and `h - G x in K`
- positive-orthant cone
- second-order/Lorentz cones
- exponential cones
- dense matrix input
- native CSC/CSR sparse storage without automatic densification on successful sparse solves
- sparse linear/SOC/exponential cone Jacobian and curvature assembly
- full sparse symmetric saddle-point Newton/KKT system
- sparse symbolic/numeric `LDL^T` factorization
- native approximate-minimum-degree (AMD-style) ordering
- reverse Cuthill-McKee fallback ordering
- cone-preserving iterative equilibration
- dynamic quasi-definite regularization
- residual-based iterative refinement
- best-iterate restoration
- ECOS inaccurate-optimal status offset
- sparse primal-unbounded and primal-infeasible ray certificates
- persistent `setup -> solve -> update -> solve -> cleanup` workflow
- symbolic-factorization and warm-start reuse across separate workspace solves
- same-pattern sparse matrix value updates without repeating symbolic analysis
- boolean/integer branch-and-bound without densifying sparse node problems
- fixed sparse bound-row structure and symbolic/warm reuse across ECOS_BB nodes
- ECOS-compatible exit/status constants and result fields
- primal/dual/slack vectors, objectives, residuals, gaps, certificates, and sparse profiling diagnostics

## Sparse Newton architecture

For a CSC problem the continuous sparse solver forms the symmetric saddle-point system directly:

```text
[ H   A^T   J^T ] [dx]     [rhs_x]
[ A   -R      0 ] [dy]  =  [rhs_y]
[ J     0   -S/L] [dl]     [rhs_l]
```

`J` is the sparse scalar-cone Jacobian and `H` contains the sparse SOC/exponential-cone curvature terms. The KKT structure is ordered and symbolically analyzed once, then numerically refactorized as cone scaling changes. Predictor and corrector solves share each numeric factorization.

With a persistent `ecos_workspace`, the symbolic elimination tree/permutation is exported to the workspace cache. A later solve reuses it when the KKT pattern matches. The previous interior point is also retained as a warm start. A change to `c`, `h`, or `b` preserves both caches; a sparse matrix value change with the same pattern preserves the symbolic cache but invalidates the warm point; a structural matrix change invalidates both as needed.

## Ordering

The default sparse ordering is now a native **approximate minimum-degree** implementation. It maintains sparse dynamic adjacency lists, uses lazy degree-heap updates, performs exact clique fill updates for ordinary fronts, and switches to a bounded-cost degree approximation for very large fronts. RCM remains selectable for comparison.

This is an AMD-style implementation, not a line-for-line port of SuiteSparse AMD's supervariable/aggressive-absorption machinery. On the regression star graph it reduces symbolic off-diagonal `L` entries from 300 with identity ordering to 24.

## Cone-preserving equilibration

Sparse problems are equilibrated by default. Variable columns are scaled iteratively. Equality rows may scale independently; positive-orthant rows may scale independently; every SOC block and every exponential-cone triple receives one shared positive row scale so cone membership is unchanged.

The final primal/dual/slack solution is mapped back to the user's original coordinates. Scale extrema are reported in `ecos_result`.

## Sparse certificates

When the main sparse solve does not converge, v0.4 first tries sparse homogeneous ray feasibility problems instead of immediately converting to dense form.

For primal unboundedness it searches for a ray `d` satisfying

```text
A*d = 0
-G*d in K
c^T*d = -1
```

For primal infeasibility it searches for a dual ray `(y,z)` satisfying

```text
A^T*y + G^T*z = 0
b^T*y + h^T*z = -1
z in K*
```

Linear and SOC blocks are self-dual. For an exponential dual block `(u,v,w)`, the implementation imposes `(-v,-u,e*w) in K_exp`. Valid normalized rays are returned in `primal_certificate` or `dual_certificate`.

This provides sparse conic certificates, but it is **not yet the integrated homogeneous self-dual embedding used internally by upstream ECOS**.

## Sparse ECOS_BB

Sparse mixed-integer problems no longer densify before branch-and-bound. Each integer variable receives a fixed pair of sparse lower/upper bound rows. Inactive rows retain explicit structural zeros, so branching changes numerical values without changing the node KKT sparsity pattern. Child nodes reuse the symbolic analysis and previous interior point. The final fixed-integer resolve also stays sparse.

## 10,000-variable structural benchmark

`example/sparse_large_lp.f90` solves a diagonal 10,000-variable LP. In the strict `-O0 -fcheck=all` validation build it reports approximately:

```text
max |x-1|:             6.25e-10
KKT upper nnz:         30000
LDL off-diagonal nnz:  10000
symbolic analyses:         1
numeric factorizations:    4
LDL/KKT fill ratio:     0.333
```

On the validation container the v0.4 executable used about **11 MB peak RSS** and **0.30 s wall time**. The same strict v0.2 example took about **0.67 s** on the same container. These are environment-specific structural checks, not universal performance claims.

## MatrixExtra integration

`MatrixExtra-fortran` is **not** a mandatory dependency of the solver. Its sparse construction/conversion API is useful, but the translated MatrixExtra direct sparse-solvers currently densify internally and therefore are not used for the ECOS KKT factorization.

An optional zero-densification adapter is provided at:

```text
integration/matrixextra-adapter/
```

It accepts Matrix/MatrixExtra CSC and COO objects and passes sparse arrays directly into ECOS.

## Remaining ECOS parity work

v0.4 is substantially closer to ECOS's sparse execution profile, but it is still not a line-for-line modern-Fortran port of the complete ECOS core. The principal remaining differences are:

- integrated homogeneous self-dual embedding rather than post-failure sparse ray problems;
- exact ECOS/Nesterov-Todd SOC scaling and ECOS exponential-cone scaling layout;
- SuiteSparse AMD's full supervariable/aggressive-absorption implementation rather than the native AMD-style ordering;
- ECOS's exact static/dynamic regularization constants and all numerical-recovery paths;
- complete best-iterate and inaccurate infeasibility status semantics;
- original ECOS_BB heuristics/cutoff logic and all node-selection rules;
- zero-allocation numeric KKT value updates; v0.4 reuses symbolic structure but still rebuilds sparse numeric KKT values each Newton iteration.

See `TRANSLATION_COVERAGE.md` for details.

## Build

```text
fpm build
fpm test
```

The core package requires only a Fortran 2018 compiler.

## Licensing

The supplied ECOSolveR package declares GPL >= 3. The referenced ECOS source and retained notices are documented in `UPSTREAM_PROVENANCE.md`. The sparse LDL implementation follows the published SuiteSparse LDL algorithm; its attribution and LGPL-2.1-or-later notice are retained in `LICENSES/` and `NOTICE.md`. The combined root package is distributed under GPL-3.0-or-later.
