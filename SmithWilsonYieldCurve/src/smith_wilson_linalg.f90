! SmithWilsonYieldCurve modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only

module smith_wilson_linalg
   use smith_wilson_kinds, only : dp
   implicit none
   private

   public :: solve_linear_system

contains

   subroutine solve_linear_system(a, b, x, info, relative_tolerance)
      real(dp), intent(in) :: a(:, :)
      real(dp), intent(in) :: b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: relative_tolerance

      real(dp), allocatable :: lu(:, :), rhs(:), work_row(:)
      real(dp) :: factor, pivot_value, scale, tolerance, tmp
      integer :: i, k, n, pivot

      info = 0
      n = size(a, 1)

      if (n < 1 .or. size(a, 2) /= n .or. size(b) /= n .or. size(x) /= n) then
         info = -1
         if (size(x) > 0) x = 0.0_dp
         return
      end if

      allocate(lu(n, n), rhs(n), work_row(n))
      lu = a
      rhs = b
      x = 0.0_dp

      scale = max(1.0_dp, maxval(sum(abs(lu), dim=2)))
      tolerance = 100.0_dp * epsilon(1.0_dp) * scale
      if (present(relative_tolerance)) then
         tolerance = max(0.0_dp, relative_tolerance) * scale
      end if

      do k = 1, n - 1
         pivot = k - 1 + maxloc(abs(lu(k:n, k)), dim=1)
         pivot_value = abs(lu(pivot, k))
         if (pivot_value <= tolerance) then
            info = k
            return
         end if

         if (pivot /= k) then
            work_row = lu(k, :)
            lu(k, :) = lu(pivot, :)
            lu(pivot, :) = work_row
            tmp = rhs(k)
            rhs(k) = rhs(pivot)
            rhs(pivot) = tmp
         end if

         do i = k + 1, n
            factor = lu(i, k) / lu(k, k)
            lu(i, k) = 0.0_dp
            lu(i, k + 1:n) = lu(i, k + 1:n) - factor * lu(k, k + 1:n)
            rhs(i) = rhs(i) - factor * rhs(k)
         end do
      end do

      if (abs(lu(n, n)) <= tolerance) then
         info = n
         return
      end if

      x(n) = rhs(n) / lu(n, n)
      do i = n - 1, 1, -1
         x(i) = (rhs(i) - dot_product(lu(i, i + 1:n), x(i + 1:n))) / lu(i, i)
      end do
   end subroutine solve_linear_system

end module smith_wilson_linalg
