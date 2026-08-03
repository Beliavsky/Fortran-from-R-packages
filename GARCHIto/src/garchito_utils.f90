! SPDX-License-Identifier: GPL-3.0-only
module garchito_utils
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchito_kinds, only : dp
   implicit none
   private

   public :: mean_value, median_value, ols_line, all_finite

contains

   pure logical function all_finite(x)
      real(dp), intent(in) :: x(:)
      all_finite = all(ieee_is_finite(x))
   end function all_finite

   pure real(dp) function mean_value(x)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         mean_value = 0.0_dp
      else
         mean_value = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   function median_value(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      real(dp), allocatable :: work(:)
      real(dp) :: key
      integer :: i, j, n

      n = size(x)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      allocate(work(n))
      work = x
      do i = 2, n
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      if (mod(n, 2) == 1) then
         value = work((n + 1) / 2)
      else
         value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
      end if
   end function median_value

   subroutine ols_line(x, y, slope, intercept, sigma)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: slope, intercept, sigma
      real(dp) :: mx, my, sxx, sxy, rss
      integer :: n

      n = size(x)
      mx = mean_value(x)
      my = mean_value(y)
      sxx = sum((x - mx)**2)
      sxy = sum((x - mx) * (y - my))
      if (sxx > epsilon(1.0_dp) * max(1.0_dp, sum(x*x))) then
         slope = sxy / sxx
      else
         slope = 0.0_dp
      end if
      intercept = my - slope * mx
      rss = sum((y - intercept - slope*x)**2)
      if (n > 2) then
         sigma = sqrt(max(rss / real(n - 2, dp), tiny(1.0_dp)))
      else if (n > 0) then
         sigma = sqrt(max(rss / real(n, dp), tiny(1.0_dp)))
      else
         sigma = 1.0_dp
      end if
      sigma = max(sigma, 1.0e-10_dp)
   end subroutine ols_line

end module garchito_utils
