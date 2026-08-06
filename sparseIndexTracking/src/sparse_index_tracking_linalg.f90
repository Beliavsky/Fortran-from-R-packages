! sparseIndexTracking modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only

module sparse_index_tracking_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use sparse_index_tracking_kinds, only : dp
   implicit none
   private

   public :: largest_eigenvalue_psd
   public :: all_finite

   interface all_finite
      module procedure all_finite_vector
      module procedure all_finite_matrix
   end interface all_finite

contains

   subroutine largest_eigenvalue_psd(a, eigenvalue, info, tolerance, max_iterations)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(out) :: eigenvalue
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_iterations

      real(dp), allocatable :: v(:), y(:)
      real(dp) :: lambda, lambda_old, norm_y, tol
      integer :: i, iteration, n, maxit

      n = size(a, 1)
      eigenvalue = 0.0_dp
      info = 0

      if (n < 1 .or. size(a, 2) /= n) then
         info = 1
         return
      end if
      if (.not. all_finite_matrix(a)) then
         info = 2
         return
      end if

      tol = 1.0e-12_dp
      if (present(tolerance)) tol = max(tolerance, epsilon(1.0_dp))
      maxit = 500
      if (present(max_iterations)) maxit = max(1, max_iterations)

      allocate(v(n), y(n))
      do i = 1, n
         v(i) = sin(real(i, dp) * 1.4142135623730950488_dp) + &
                cos(real(i, dp) * 0.6180339887498948482_dp)
      end do
      norm_y = norm2(v)
      if (norm_y <= tiny(1.0_dp)) then
         v = 1.0_dp / sqrt(real(n, dp))
      else
         v = v / norm_y
      end if

      lambda_old = -huge(1.0_dp)
      do iteration = 1, maxit
         y = matmul(a, v)
         norm_y = norm2(y)
         if (norm_y <= 100.0_dp * tiny(1.0_dp)) then
            eigenvalue = 0.0_dp
            return
         end if
         v = y / norm_y
         lambda = dot_product(v, matmul(a, v))
         if (iteration > 1) then
            if (abs(lambda - lambda_old) <= tol * max(1.0_dp, abs(lambda))) then
               eigenvalue = max(lambda, 0.0_dp)
               return
            end if
         end if
         lambda_old = lambda
      end do

      eigenvalue = max(lambda_old, 0.0_dp)
      info = 3
   end subroutine largest_eigenvalue_psd


   pure function all_finite_vector(x) result(ok)
      real(dp), intent(in) :: x(:)
      logical :: ok

      ok = all(ieee_is_finite(x))
   end function all_finite_vector


   pure function all_finite_matrix(x) result(ok)
      real(dp), intent(in) :: x(:, :)
      logical :: ok

      ok = all(ieee_is_finite(x))
   end function all_finite_matrix

end module sparse_index_tracking_linalg
