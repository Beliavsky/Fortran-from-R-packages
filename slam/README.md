# slam-fortran

A modern Fortran/FPM translation of the numerical computational core of the R
package **slam 0.1-56** (Sparse Lightweight Arrays and Matrices), by Kurt
Hornik, David Meyer, and Christian Buchta.

The original package is licensed under **GPL-2**. This translation keeps the
same license family; see `COPYING`, `LICENSE.note`, and `ORIGINAL_DESCRIPTION`.
Every translated source file carries an SPDX `GPL-2.0-only` header.

## Scope

The port translates the algorithms and data structures that are useful in a
statically typed numerical Fortran program. It does not try to reproduce R's
S3/S4 dispatch machinery or R-specific object metadata.

Implemented computational areas:

- `simple_triplet_matrix` sparse matrix representation with 1-based indices
- validation and duplicate-index checking
- dense <-> triplet conversion
- zero and diagonal constructors
- transpose and dimension reshaping in column-major order
- row/column sums and means, with optional NaN removal
- row/column p-norms, including infinity norm
- matrix mean
- sparse addition, subtraction, Hadamard multiplication, scalar
  multiplication/division, and positive scalar powers
- row and column scaling, including non-finite-value handling
- row/column extraction and coordinate replacement with value recycling
- row/column binding
- sparse-sparse `tcrossprod`, `crossprod`, and matrix product
- sparse-dense and dense-sparse products
- `simple_sparse_array` representation with arbitrary rank
- dense-flat <-> sparse-array conversion using Fortran/R column-major order
- sparse-array validation, permutation, reshape, reduction, dimension dropping,
  extension, binding, assignment, and matrix conversion
- grouped sum rollups for triplet matrices and sparse arrays
- column/row/cross/tcross apply operations with scalar-return Fortran callbacks
- CLUTO sparse matrix file input/output
- MC sparse matrix file input/output
- reusable column-major index conversion and sorting utilities corresponding to
  the low-level C indexing/grouping machinery

There is no substantive plotting code in `slam`; nothing graphical was ported.

## Deliberate interface differences

The R package supports dynamically typed payloads and extensive R object
semantics. This Fortran port uses `real(dp)` numerical payloads. In particular,
the following are intentionally outside the translated numerical API:

- S3/S4 method dispatch, printing, and R option plumbing
- `dimnames`, names, factors, and other R attributes
- coercion adapters for Matrix, SparseM, and spam classes
- character, raw, expression, list, and general complex-valued sparse payloads
- arbitrary-length return values from R callbacks; the Fortran apply API uses
  scalar-return procedures
- R's full logical/character/negative/missing subscript grammar; the Fortran API
  exposes explicit integer row/column or coordinate operations
- the most dynamic `rollup(..., FUN=...)` expansion modes; the optimized sum
  path is translated directly, while arbitrary reductions can be expressed with
  Fortran callback code around the sparse containers
- `simplify_simple_sparse_array`, whose purpose is to reshape list-valued R
  cells and therefore has no direct `real(dp)` analogue

NaN is used for R-like numeric missing values where applicable. The product
routines detect non-finite stored values and use a dense fallback when needed so
that cases such as `0 * Inf` retain IEEE semantics instead of being silently
lost by sparse traversal.

## Source mapping

| Original source | Fortran translation |
| --- | --- |
| `R/matrix.R`, `R/stm.R`, `src/sparse.c` | `src/slam_stm.f90` |
| `R/array.R`, `R/reduce.R`, `R/subassign.R` | `src/slam_ssa.f90` |
| `R/abind.R` | `src/slam_ssa.f90` |
| `R/crossprod.R`, cross-product C kernels | `src/slam_stm.f90` |
| `R/apply.R`, `src/apply.c` | `src/slam_apply.f90` |
| `R/rollup.R`, `src/grouped.c` | `src/slam_rollup.f90` |
| `R/foreign.R` | `src/slam_io.f90` |
| indexing/matching portions of `src/util.c` | `src/slam_utils.f90` |
| R namespace/class glue and `slam_options` | intentionally omitted |

## Building with FPM

```text
fpm build
fpm test
fpm run --example basic
```

The included `fpm.toml` uses free-form modern Fortran and disables implicit
typing and implicit external procedures.

FPM was not installed in the translation environment, so validation here used
the equivalent direct GNU Fortran build with GNU Fortran 14.2.0:

```text
gfortran -std=f2018 -Wall -Wextra -Wconversion -fcheck=all -fbacktrace -O0 \
  -c src/slam_kinds.f90 src/slam_utils.f90 src/slam_stm.f90 \
     src/slam_ssa.f90 src/slam_rollup.f90 src/slam_apply.f90 \
     src/slam_io.f90 src/slam.f90

gfortran -std=f2018 -fcheck=all test/test_slam.f90 *.o -o test_slam
./test_slam
```

The bounds-checked tests pass. Compiler warnings about exact real comparisons
are intentional: sparse storage must distinguish exact zero from nonzero and
NaN values.

## Minimal example

```fortran
program basic
    use slam
    implicit none
    real(dp) :: a(3,3)
    real(dp), allocatable :: gram(:,:)
    type(simple_triplet_matrix) :: x

    a = reshape([1.0_dp,0.0_dp,2.0_dp, &
                 0.0_dp,3.0_dp,0.0_dp, &
                 4.0_dp,0.0_dp,5.0_dp], [3,3])

    x = dense_to_stm(a)
    gram = stm_crossprod(x)

    print '(a,i0)', 'nonzeros: ', x%nnz()
    print '(3f10.3)', gram
end program basic
```

## Representation notes

Both sparse representations use 1-based indices, matching both R and normal
Fortran indexing. Dense linearization uses column-major ordering. A sparse array
is represented by an `(nnz, rank)` coordinate matrix, a `real(dp)` value vector,
and its dimension vector.
