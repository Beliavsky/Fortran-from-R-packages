! SPDX-License-Identifier: GPL-3.0-only
module matrix_sparse_stats
   use matrix_kinds, only : dp
   use matrix_sparse, only : csr_matrix, csr_transpose, csr_multiply
   use matrix_status, only : matrix_success
   implicit none
   private
   public :: csr_row_sums, csr_col_sums, csr_row_means, csr_col_means
   public :: csr_trace, csr_frobenius_norm, csr_one_norm, csr_infinity_norm
   public :: csr_crossprod, csr_tcrossprod

contains

   function csr_row_sums(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: x(:)
      integer :: i
      allocate(x(a%nrow), source=0.0_dp)
      do i = 1, a%nrow
         x(i) = sum(a%values(a%row_ptr(i):a%row_ptr(i + 1) - 1))
      end do
   end function csr_row_sums

   function csr_col_sums(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: x(:)
      integer :: k
      allocate(x(a%ncol), source=0.0_dp)
      do k = 1, a%nnz()
         x(a%col_ind(k)) = x(a%col_ind(k)) + a%values(k)
      end do
   end function csr_col_sums

   function csr_row_means(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: x(:)
      x = csr_row_sums(a)
      if (a%ncol > 0) x = x / real(a%ncol, dp)
   end function csr_row_means

   function csr_col_means(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: x(:)
      x = csr_col_sums(a)
      if (a%nrow > 0) x = x / real(a%nrow, dp)
   end function csr_col_means

   real(dp) function csr_trace(a) result(x)
      type(csr_matrix), intent(in) :: a
      integer :: i, k
      x = 0.0_dp
      do i = 1, min(a%nrow, a%ncol)
         do k = a%row_ptr(i), a%row_ptr(i + 1) - 1
            if (a%col_ind(k) == i) then
               x = x + a%values(k)
               exit
            end if
         end do
      end do
   end function csr_trace

   real(dp) function csr_frobenius_norm(a) result(x)
      type(csr_matrix), intent(in) :: a
      x = sqrt(sum(a%values * a%values))
   end function csr_frobenius_norm

   real(dp) function csr_one_norm(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp), allocatable :: sums(:)
      integer :: k
      allocate(sums(a%ncol), source=0.0_dp)
      do k = 1, a%nnz()
         sums(a%col_ind(k)) = sums(a%col_ind(k)) + abs(a%values(k))
      end do
      if (a%ncol == 0) then
         x = 0.0_dp
      else
         x = maxval(sums)
      end if
   end function csr_one_norm

   real(dp) function csr_infinity_norm(a) result(x)
      type(csr_matrix), intent(in) :: a
      real(dp) :: row_sum
      integer :: i
      x = 0.0_dp
      do i = 1, a%nrow
         row_sum = sum(abs(a%values(a%row_ptr(i):a%row_ptr(i + 1) - 1)))
         x = max(x, row_sum)
      end do
   end function csr_infinity_norm

   subroutine csr_crossprod(a, c, info)
      type(csr_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: c
      integer, intent(out) :: info
      type(csr_matrix) :: at
      call csr_transpose(a, at)
      call csr_multiply(at, a, c, info)
   end subroutine csr_crossprod

   subroutine csr_tcrossprod(a, c, info)
      type(csr_matrix), intent(in) :: a
      type(csr_matrix), intent(out) :: c
      integer, intent(out) :: info
      type(csr_matrix) :: at
      call csr_transpose(a, at)
      call csr_multiply(a, at, c, info)
   end subroutine csr_tcrossprod

end module matrix_sparse_stats
