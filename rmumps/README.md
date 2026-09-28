# rmumps - Fortran translation

This package translates the computational interface of the R package `rmumps` to modern free-form Fortran with FPM packaging.

The current implementation is portable and has no external numerical-library dependency. It stores the input sparse matrix as triplets, assembles a dense matrix, and caches an internal LU factorization with partial pivoting. The translated API supports dense and sparse-triplet right-hand sides, transposed solves, inverse and inverse-transpose requests, determinant calculations, triplet access, value updates, symmetry modes, and permutation-code state.

This is **not yet a replacement for MUMPS on large sparse systems**. The current backend is dense and serial. It exposes MUMPS-shaped control/information state and reproduces selected observable control semantics, but it does not implement MUMPS ordering algorithms, multifrontal sparse factorization, most MUMPS-specific control effects and diagnostics, distributed-memory behavior, or Rcpp object semantics.

## Build

```text
fpm build
fpm test
```

A direct compiler check can also be run with a recent gfortran:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all -J build -I build \
  src/rmumps_kinds.f90 src/rmumps_api.f90 test/test_rmumps.f90 -o build/test_rmumps
```

## Main API

The public module is `rmumps_api`. Important entry points include:

- `rmumps_init_triplet` - create a matrix from one-based sparse triplets.
- `rmumps_init_triplet_zero_based` - create a matrix from upstream/Rcpp-style zero-based triplets.
- `rmumps_symbolic`, `rmumps_numeric` - explicitly mark symbolic analysis and compute/cache the numeric LU factorization.
- `rmumps_solve` - generic solve for a dense vector, dense matrix, or sparse-triplet RHS.
- `rmumps_solve_inplace` - generic in-place solve for dense vectors and matrices.
- `rmumps_solvet` - generic solve with `transpose(A)` using the same cached LU factors.
- `rmumps_inverse`, `rmumps_inverse_transpose` - dense inverse operations.
- `rmumps_det`, `rmumps_log2det`, `rmumps_determinant` - determinant interfaces.
- `rmumps_set_values` - replace numerical values while preserving the triplet sparsity pattern.
- `rmumps_get_triplet` - retrieve cloned triplet arrays.
- `rmumps_set_permutation`, `rmumps_get_permutation` - store/query rmumps ordering codes.
- `rmumps_set_perm_in` - validate and store a user-supplied symmetric permutation request.
- `rmumps_set_copy`, `rmumps_get_copy` - retain the R-wrapper copy flag as compatibility state.
- `rmumps_set_ncore`, `rmumps_get_ncore` - retain the requested core count; the fallback remains serial.
- `rmumps_set_rhs`, `rmumps_get_rhs`, `rmumps_set_mrhs`, `rmumps_get_mrhs` - maintain public RHS state.
- `rmumps_set_icntl`, `rmumps_get_icntl` - maintain the 33-entry MUMPS-style integer-control state.
- `rmumps_set_cntl`, `rmumps_get_cntl` - maintain the 5-entry MUMPS-style real-control state.
- `rmumps_set_keep`, `rmumps_get_keep` - maintain the 500-entry MUMPS-style KEEP state.
- `rmumps_get_infos` - expose MUMPS-shaped INFO/RINFO/INFOG/RINFOG arrays with fallback status and determinant data.
- `rmumps_dim`, `rmumps_nrow`, `rmumps_ncol` - matrix dimensions.

Triplet indices are one-based in the primary Fortran API. `rmumps_init_triplet_zero_based` is provided for direct parity with the upstream Rcpp constructor. Duplicate triplet entries are summed during dense assembly.

For `sym=1` or `sym=2`, off-diagonal entries are mirrored to the opposite triangle. As with upstream rmumps, callers should supply only one triangle of a symmetric matrix rather than both triangles.

Sparse RHS routines accept one-based `(i,j,v)` triplets and sum duplicate entries. Passing an empty `0 x 0` sparse RHS requests the inverse, matching the special case used by upstream rmumps; the analogous `rmumps_solvet` call returns the inverse transpose.

## Example

```fortran
use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_solve

type(rmumps_matrix) :: a
integer :: info
integer :: i(4), j(4)
real(dp) :: v(4), b(2), x(2)

i = [1, 1, 2, 2]
j = [1, 2, 1, 2]
v = [3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp]
b = [5.0_dp, 5.0_dp]

call rmumps_init_triplet(a, i, j, v, 2, info=info)
call rmumps_solve(a, b, x, info)
```

## Validation

Deterministic tests cover:

- vector and multiple dense right-hand sides;
- in-place solves;
- transposed solves that reuse cached LU factors, including a matrix requiring pivot swaps;
- sparse-triplet RHS matrices and duplicate RHS entries;
- inverse and inverse-transpose requests through the upstream `0 x 0` sparse-RHS convention;
- determinant sign and logarithmic/non-logarithmic modulus;
- symmetric matrix input;
- duplicate matrix triplets;
- MUMPS-compatible `-10` status for singular matrices and input validation;
- triplet retrieval, value replacement, dimensions, and permutation-code state;
- upstream rmumps control defaults, including tested KEEP integer-size entries, plus ICNTL/CNTL/KEEP set/get behavior;
- observable `ICNTL(9)` transpose-solve and `ICNTL(20)` dense/sparse RHS-mode semantics;
- MUMPS-shaped information arrays and determinant mantissa/exponent encoding;
- copy/core-count compatibility state, stored RHS state, and user-permutation validation/invalidation.

The source and tests are checked with strict gfortran warnings and runtime bounds/checking flags. FPM is the intended package build system.

## Translation coverage

Package status: **partial**.

Coverage: **13 of 13 (100%)** exported computational R functions are mapped. This fraction measures mapped computational functions under the requested counting convention; it does **not** mean complete R compatibility or complete MUMPS feature compatibility. In particular, the sparse multifrontal MUMPS backend remains untranslated/replaced by the current dense fallback.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `Rmumps` | `rmumps_api` | `rmumps_init_triplet`, `rmumps_init_triplet_zero_based` | partial |
| `Rmumps__solveptr` | `rmumps_api` | `rmumps_solve_vector_inplace`, `rmumps_solve_matrix_inplace` | substantial |
| `Rmumps__ptr_ijv` | `rmumps_api` | `rmumps_init_triplet_zero_based`, `rmumps_init_triplet` | substantial |
| `Rmumps__triplet` | `rmumps_api` | `rmumps_get_triplet` | substantial |
| `Rmumps__set_mat_ptr` | `rmumps_api` | `rmumps_set_values` | substantial |
| `Rmumps__set_permutation` | `rmumps_api` | `rmumps_set_permutation` | partial |
| `Rmumps__get_permutation` | `rmumps_api` | `rmumps_get_permutation` | complete |
| `solve.Rcpp_Rmumps` | `rmumps_api` | `rmumps_solve_vector`, `rmumps_solve_matrix`, `rmumps_solve_sparse_rhs`, `rmumps_inverse` | substantial |
| `solvet.Rcpp_Rmumps` | `rmumps_api` | `rmumps_solvet_vector`, `rmumps_solvet_matrix`, `rmumps_solvet_sparse_rhs`, `rmumps_inverse_transpose` | substantial |
| `determinant.Rcpp_Rmumps` | `rmumps_api` | `rmumps_determinant` | substantial |
| `dim.Rcpp_Rmumps` | `rmumps_api` | `rmumps_dim` | complete |
| `nrow.Rcpp_Rmumps` | `rmumps_api` | `rmumps_nrow` | complete |
| `ncol.Rcpp_Rmumps` | `rmumps_api` | `rmumps_ncol` | complete |

See `API_COVERAGE.md` for the counting basis and class-method compatibility notes, and `NOTICE.md` for provenance and backend limitations.
