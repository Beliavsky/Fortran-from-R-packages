# Changelog

## 0.4.0

- Added native approximate-minimum-degree (AMD-style) sparse ordering; RCM remains available.
- Added cone-preserving iterative sparse equilibration with variable/equality/cone scale tracking.
- Added dynamic KKT regularization retry schedule.
- Added residual-improvement iterative refinement profiling.
- Added best-iterate restoration and `ECOS_OPTIMAL + ECOS_INACC_OFFSET` handling.
- Added persistent sparse symbolic-factor cache in `ecos_workspace`.
- Added persistent interior-point warm starts across workspace solves.
- `c`, `h`, and `b` updates preserve symbolic and warm caches.
- Same-pattern sparse matrix value updates preserve symbolic analysis and invalidate only the warm point.
- Structural sparse matrix updates invalidate symbolic analysis.
- Added sparse homogeneous primal-unbounded ray certificates.
- Added sparse dual-ray primal-infeasibility certificates.
- Added the dual exponential-cone transform used by certificate construction.
- Replaced the large sparse diagnostic densification path with sparse certificates; optional dense fallback is limited to small problems.
- Sparse ECOS_BB no longer densifies node problems.
- Added fixed structural branch-bound rows and symbolic/warm reuse across nodes.
- Final fixed-integer sparse re-solve stays sparse.
- Added sparse factor/order/refinement timing and fill-ratio diagnostics.
- Added scale-range, cache-reuse, regularization-update, and certificate fields to `ecos_result`.
- Extended `ecos_control` with AMD, equilibration, regularization, and certificate controls.
- Added `test_v04_features` covering AMD fill, equilibration, persistent cache reuse, sparse certificates, exponential dual certificates, sparse ECOS_BB reuse, and inaccurate-optimal status.
- Added `workspace_reuse` example.
- Updated the MatrixExtra adapter project to version 0.4.0.

## 0.2.0

- CSC input remains sparse in continuous solves.
- Added native CSR storage and CSC/CSR sparse matvec operations.
- Added sparse linear/SOC/exponential cone Jacobian and curvature assembly.
- Added sparse symmetric saddle-point KKT assembly.
- Added SuiteSparse-LDL-style symbolic/numeric `LDL^T` factorization.
- Added RCM ordering and iterative refinement.
- Reused symbolic analysis over IPM iterations and numeric factors for predictor/corrector solves.
- Added sparse diagnostics and 2,000/10,000-variable sparse regressions/examples.
- Added optional MatrixExtra adapter.

## 0.1.0

- Initial modern Fortran/FPM translation.
- Positive orthant, SOC, and exponential cone support.
- Dense and CSC problem input.
- Predictor/corrector primal-dual solver.
- Mixed-integer/boolean branch-and-bound.
- Persistent setup/update/solve lifecycle.
- ECOS-compatible settings, status constants, and result structure.
