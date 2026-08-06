! sparseIndexTracking modern Fortran translation
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-only

module sparse_index_tracking_projection
   use sparse_index_tracking_kinds, only : dp
   use sparse_index_tracking_linalg, only : all_finite
   implicit none
   private

   public :: project_capped_simplex
   public :: bisection

contains

   subroutine project_capped_simplex(z, upper_bound, weights, info)
      real(dp), intent(in) :: z(:)
      real(dp), intent(in) :: upper_bound
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out) :: info

      real(dp), allocatable :: c(:)

      allocate(c(size(z)))
      c = -2.0_dp * z
      call bisection(c, upper_bound, weights, info)
   end subroutine project_capped_simplex


   subroutine bisection(c, upper_bound, weights, info)
      real(dp), intent(in) :: c(:)
      real(dp), intent(in) :: upper_bound
      real(dp), allocatable, intent(out) :: weights(:)
      integer, intent(out) :: info

      real(dp) :: difference, mu, mu_high, mu_low, total
      integer :: i, iteration, n

      n = size(c)
      allocate(weights(n))
      weights = 0.0_dp
      info = 0

      if (n < 1 .or. upper_bound <= 0.0_dp .or. .not. all_finite(c)) then
         info = 1
         return
      end if
      if (real(n, dp) * upper_bound < 1.0_dp - 100.0_dp * epsilon(1.0_dp)) then
         info = 2
         return
      end if

      mu_low = minval(-c - 2.0_dp * upper_bound) - 1.0_dp
      mu_high = maxval(-c) + 1.0_dp

      do iteration = 1, 200
         mu = 0.5_dp * (mu_low + mu_high)
         weights = min(upper_bound, max(0.0_dp, -0.5_dp * (mu + c)))
         total = sum(weights)
         if (total > 1.0_dp) then
            mu_low = mu
         else
            mu_high = mu
         end if
         if (abs(total - 1.0_dp) <= 20.0_dp * epsilon(1.0_dp)) exit
      end do

      difference = 1.0_dp - sum(weights)
      if (difference > 0.0_dp) then
         do i = 1, n
            if (weights(i) < upper_bound) then
               total = min(difference, upper_bound - weights(i))
               weights(i) = weights(i) + total
               difference = difference - total
               if (difference <= 50.0_dp * epsilon(1.0_dp)) exit
            end if
         end do
      else if (difference < 0.0_dp) then
         do i = 1, n
            if (weights(i) > 0.0_dp) then
               total = min(-difference, weights(i))
               weights(i) = weights(i) - total
               difference = difference + total
               if (-difference <= 50.0_dp * epsilon(1.0_dp)) exit
            end if
         end do
      end if

      if (abs(sum(weights) - 1.0_dp) > 1.0e-11_dp .or. &
          any(weights < -1.0e-13_dp) .or. &
          any(weights > upper_bound + 1.0e-13_dp)) info = 3
   end subroutine bisection

end module sparse_index_tracking_projection
