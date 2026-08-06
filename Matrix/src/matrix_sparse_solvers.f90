! SPDX-License-Identifier: GPL-3.0-only
module matrix_sparse_solvers
   use matrix_kinds, only : dp
   use matrix_status, only : matrix_success, matrix_err_shape
   use matrix_sparse, only : csr_matrix, csr_to_dense
   use matrix_decompositions, only : lu_factor, lu_solve, cholesky_factor, cholesky_solve, least_squares
   implicit none
   private

   type, public :: sparse_lu_factor
      integer :: n = 0
      real(dp), allocatable :: lu(:,:)
      integer, allocatable :: piv(:)
   end type sparse_lu_factor

   type, public :: sparse_cholesky_factor
      integer :: n = 0
      real(dp), allocatable :: l(:,:)
   end type sparse_cholesky_factor

   public :: csr_lu_factorize, csr_lu_solve
   public :: csr_cholesky_factorize, csr_cholesky_solve
   public :: csr_least_squares

contains

   subroutine csr_lu_factorize(a, factor, info)
      type(csr_matrix), intent(in) :: a
      type(sparse_lu_factor), intent(out) :: factor
      integer, intent(out) :: info
      real(dp), allocatable :: dense(:,:)
      if (a%nrow /= a%ncol) then
         factor%n = 0
         allocate(factor%lu(0, 0), factor%piv(0))
         info = matrix_err_shape
         return
      end if
      dense = csr_to_dense(a)
      call lu_factor(dense, factor%lu, factor%piv, info)
      if (info == matrix_success) factor%n = a%nrow
   end subroutine csr_lu_factorize

   subroutine csr_lu_solve(factor, b, x, info)
      type(sparse_lu_factor), intent(in) :: factor
      real(dp), intent(in) :: b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      call lu_solve(factor%lu, factor%piv, b, x, info)
   end subroutine csr_lu_solve

   subroutine csr_cholesky_factorize(a, factor, info, tol)
      type(csr_matrix), intent(in) :: a
      type(sparse_cholesky_factor), intent(out) :: factor
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: dense(:,:)
      if (a%nrow /= a%ncol) then
         factor%n = 0
         allocate(factor%l(0, 0))
         info = matrix_err_shape
         return
      end if
      dense = csr_to_dense(a)
      call cholesky_factor(dense, factor%l, info, tol)
      if (info == matrix_success) factor%n = a%nrow
   end subroutine csr_cholesky_factorize

   subroutine csr_cholesky_solve(factor, b, x, info)
      type(sparse_cholesky_factor), intent(in) :: factor
      real(dp), intent(in) :: b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      call cholesky_solve(factor%l, b, x, info)
   end subroutine csr_cholesky_solve

   subroutine csr_least_squares(a, b, x, info)
      type(csr_matrix), intent(in) :: a
      real(dp), intent(in) :: b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: dense(:,:)
      dense = csr_to_dense(a)
      call least_squares(dense, b, x, info)
   end subroutine csr_least_squares

end module matrix_sparse_solvers
