! SPDX-License-Identifier: GPL-2.0-or-later
module ceoptim_linalg
   use ceoptim_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   implicit none
   private
   public :: symmetric_eigen, covariance_factor, symmetric_pinv

contains

   subroutine covariance_factor(sigma, factor, info)
      real(dp), intent(in) :: sigma(:, :)
      real(dp), allocatable, intent(out) :: factor(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: eval(:), evec(:, :)
      real(dp) :: scale, tol
      integer :: n, j

      call symmetric_eigen(sigma, eval, evec, info)
      if (info /= 0) then
         allocate(factor(0, 0))
         return
      end if
      n = size(sigma, 1)
      scale = max(1.0_dp, maxval(abs(eval)))
      tol = 100.0_dp * epsilon(1.0_dp) * scale
      if (minval(eval) < -tol) then
         info = -2
         allocate(factor(0, 0))
         return
      end if
      allocate(factor(n, n))
      do j = 1, n
         factor(:, j) = evec(:, j) * sqrt(max(eval(j), 0.0_dp))
      end do
      info = 0
   end subroutine covariance_factor

   subroutine symmetric_pinv(a, ainv, info)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: ainv(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: eval(:), evec(:, :), tmp(:, :)
      real(dp) :: scale, tol
      integer :: n, j

      call symmetric_eigen(a, eval, evec, info)
      if (info /= 0) then
         allocate(ainv(0, 0))
         return
      end if
      n = size(a, 1)
      scale = max(1.0_dp, maxval(abs(eval)))
      tol = 100.0_dp * epsilon(1.0_dp) * real(max(1, n), dp) * scale
      allocate(tmp(n, n), ainv(n, n))
      tmp = 0.0_dp
      do j = 1, n
         if (abs(eval(j)) > tol) then
            tmp(:, j) = evec(:, j) / eval(j)
         end if
      end do
      ainv = matmul(tmp, transpose(evec))
      info = 0
   end subroutine symmetric_pinv

end module ceoptim_linalg
