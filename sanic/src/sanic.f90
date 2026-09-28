! SPDX-License-Identifier: GPL-3.0-or-later
! Public API for the Fortran translation of the R package sanic.
module sanic
   use sanic_kinds, only: dp
   use sanic_sparse, only: sparse_matrix, sparsify, sparse_from_triplets, sparse_to_dense, &
                           sparse_matvec, sparse_transpose_matvec, sparse_matmat, sparse_is_symmetric
   use sanic_direct, only: solve2_dense, solve2_sparse, solve_chol_dense, solve_chol_sparse, &
                           solve_lu_dense, solve_lu_sparse, solve_qr_dense, solve_qr_sparse, &
                           dense_is_symmetric
   use sanic_iterative, only: sanic_bicgstab, sanic_lscg, sanic_cg, solve_cg_sparse
   use sanic_spectral, only: eigen_result, svd_result, sanic_svd_bdc, sanic_svd_jacobi, &
                             eigen2, svd2
   use sanic_krylov, only: krylov_result, arnoldi_sparse, lanczos_sparse
   implicit none
   private

   public :: dp
   public :: sparse_matrix, sparsify, sparse_from_triplets, sparse_to_dense
   public :: sparse_matvec, sparse_transpose_matvec, sparse_matmat, sparse_is_symmetric
   public :: solve2_dense, solve2_sparse
   public :: solve_chol_dense, solve_chol_sparse
   public :: solve_lu_dense, solve_lu_sparse
   public :: solve_qr_dense, solve_qr_sparse
   public :: dense_is_symmetric
   public :: sanic_bicgstab, sanic_lscg, sanic_cg, solve_cg_sparse
   public :: eigen_result, svd_result, sanic_svd_bdc, sanic_svd_jacobi, eigen2, svd2
   public :: krylov_result, arnoldi_sparse, lanczos_sparse

   interface solve2
      module procedure solve2_dense
      module procedure solve2_sparse
   end interface solve2

   interface solve_chol
      module procedure solve_chol_dense
      module procedure solve_chol_sparse
   end interface solve_chol

   interface solve_lu
      module procedure solve_lu_dense
      module procedure solve_lu_sparse
   end interface solve_lu

   interface solve_qr
      module procedure solve_qr_dense
      module procedure solve_qr_sparse
   end interface solve_qr

   public :: solve2, solve_chol, solve_lu, solve_qr

end module sanic
