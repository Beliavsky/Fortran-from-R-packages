# Changelog

## 0.2.0

- Replaced the v0.1.0 dense KKT factorization with native sparse QDLDL.
- Added `scs_qdldl`, translated from the bundled Apache-2.0 QDLDL source.
- Added direct sparse upper-KKT CSC assembly.
- Reuse QDLDL symbolic analysis and allocated workspaces across diagonal
  refactorizations caused by adaptive scaling.
- Added sparse linear-system statistics to `scs_info`:
  `kkt_nnz`, `factor_nnz`, `factorizations`, and `symbolic_analyses`.
- Added direct sparse-factor solve/refactorization regression tests.
- Added a sparse-vs-dense LDL benchmark program.
- Updated license/provenance documentation for translated QDLDL.

Remaining performance item: native fill-reducing symmetric ordering (upstream
uses SuiteSparse AMD).

## 0.1.0

- Initial modern Fortran/FPM translation of the computational core of R `scs`.
- Implemented SCS HSDE/Douglas-Rachford solver, cone projections,
  normalization, Anderson acceleration, CSC operations, and a portable dense
  KKT LDL^T backend.
