# Changelog

## 0.1.1

- Fixed an FPM/gfortran compatibility failure in `clue_sumt.f90` where calls to optional, host-associated gradient procedure dummies from an internal procedure could be diagnosed as having an implicit interface by some GNU Fortran versions.
- Split SUMT gradient evaluation into module procedures with explicit `procedure(grad_fun)` interfaces for analytic gradients and a separate finite-difference fallback.
- Revalidated all eight tests with `-Werror=implicit-interface`, bounds checking, and `-O2`.

## 0.1.0

- Initial modern Fortran/FPM translation of the computational core of `clue`
  0.3-68.
- Added LSAP, hard/soft partition representations, partition lattice operations,
  agreement and dissimilarity measures.
- Added DWH and alternating Euclidean/Manhattan consensus algorithms.
- Added LP-backed Mallows/CSSD, exact k-medoids, and soft-L1 consensus using the
  supplied `lpSolve-fortran` dependency.
- Translated native ultrametric/additive-tree deviation, gradient, iterative
  projection and iterative reduction kernels.
- Added SUMT, high-level least-squares/L1 tree fitting, target-topology fits,
  PAVA and sum-of-ultrametrics fitting.
- Added Euclidean hard/fuzzy prototype clustering and validity measures.
- Added bounds-checked regression suite and examples.
