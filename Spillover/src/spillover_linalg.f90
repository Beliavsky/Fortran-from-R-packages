! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from the Spillover R package.
module spillover_linalg
   use spillover_kinds, only : dp
   use spillover_status, only : spillover_success, spillover_invalid_argument, &
      spillover_singular_matrix, spillover_not_positive_definite, set_status
   implicit none
   private

   public :: solve_linear_system
   public :: cholesky_lower
   public :: identity_matrix
   public :: is_finite_matrix

contains

   subroutine solve_linear_system(a, b, x, info, message)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:, :)
      real(dp), allocatable, intent(out) :: x(:, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp), allocatable :: aug(:, :), row_tmp(:)
      real(dp) :: pivot_abs, factor, scale
      integer :: n, nrhs, i, j, k, pivot

      call set_status(info, message, spillover_success, 'success')
      n = size(a, 1)
      nrhs = size(b, 2)
      if (size(a, 2) /= n .or. size(b, 1) /= n .or. n < 1) then
         allocate(x(0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'linear-system dimensions are inconsistent')
         return
      end if

      allocate(aug(n, n + nrhs), row_tmp(n + nrhs), x(n, nrhs))
      aug(:, 1:n) = a
      aug(:, n + 1:n + nrhs) = b
      scale = max(1.0_dp, maxval(abs(a)))

      do k = 1, n
         pivot = k
         pivot_abs = abs(aug(k, k))
         do i = k + 1, n
            if (abs(aug(i, k)) > pivot_abs) then
               pivot = i
               pivot_abs = abs(aug(i, k))
            end if
         end do
         if (pivot_abs <= epsilon(1.0_dp) * scale * real(max(1, n), dp)) then
            x = 0.0_dp
            call set_status(info, message, spillover_singular_matrix, &
               'linear system is singular to working precision')
            return
         end if
         if (pivot /= k) then
            row_tmp = aug(k, :)
            aug(k, :) = aug(pivot, :)
            aug(pivot, :) = row_tmp
         end if

         do i = k + 1, n
            factor = aug(i, k) / aug(k, k)
            aug(i, k) = 0.0_dp
            do j = k + 1, n + nrhs
               aug(i, j) = aug(i, j) - factor * aug(k, j)
            end do
         end do
      end do

      do j = 1, nrhs
         do i = n, 1, -1
            x(i, j) = aug(i, n + j)
            if (i < n) x(i, j) = x(i, j) - dot_product(aug(i, i + 1:n), x(i + 1:n, j))
            x(i, j) = x(i, j) / aug(i, i)
         end do
      end do
   end subroutine solve_linear_system

   subroutine cholesky_lower(a, l, info, message)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: l(:, :)
      integer, intent(out), optional :: info
      character(len=*), intent(out), optional :: message

      real(dp) :: value, scale
      integer :: n, i, j

      call set_status(info, message, spillover_success, 'success')
      n = size(a, 1)
      if (size(a, 2) /= n .or. n < 1) then
         allocate(l(0, 0))
         call set_status(info, message, spillover_invalid_argument, &
            'Cholesky input must be a nonempty square matrix')
         return
      end if

      allocate(l(n, n))
      l = 0.0_dp
      scale = max(1.0_dp, maxval(abs(a)))
      do i = 1, n
         do j = 1, i
            value = a(i, j)
            if (j > 1) value = value - dot_product(l(i, 1:j - 1), l(j, 1:j - 1))
            if (i == j) then
               if (value <= epsilon(1.0_dp) * scale * real(n, dp)) then
                  l = 0.0_dp
                  call set_status(info, message, spillover_not_positive_definite, &
                     'matrix is not positive definite')
                  return
               end if
               l(i, j) = sqrt(value)
            else
               l(i, j) = value / l(j, j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n, n)
      integer :: i

      a = 0.0_dp
      do i = 1, n
         a(i, i) = 1.0_dp
      end do
   end function identity_matrix

   pure logical function is_finite_matrix(a)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: a(:, :)

      is_finite_matrix = all(ieee_is_finite(a))
   end function is_finite_matrix

end module spillover_linalg
