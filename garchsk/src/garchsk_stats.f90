! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use garchsk_kinds, only : dp
   implicit none
   private
   public :: mean_value, sample_variance, sample_skewness, sample_kurtosis
   public :: covariance, all_finite

contains

   pure real(dp) function mean_value(x) result(value)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x) / real(size(x), dp)
      end if
   end function mean_value

   pure real(dp) function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      if (size(x) < 2) then
         value = 0.0_dp
         return
      end if
      m = mean_value(x)
      value = sum((x - m)**2) / real(size(x) - 1, dp)
   end function sample_variance

   pure real(dp) function sample_skewness(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m, sd
      if (size(x) < 2) then
         value = 0.0_dp
         return
      end if
      m = mean_value(x)
      sd = sqrt(max(sample_variance(x), 0.0_dp))
      if (sd <= tiny(1.0_dp)) then
         value = 0.0_dp
      else
         value = sum(((x - m) / sd)**3) / real(size(x), dp)
      end if
   end function sample_skewness

   pure real(dp) function sample_kurtosis(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: m, sd
      if (size(x) < 2) then
         value = 0.0_dp
         return
      end if
      m = mean_value(x)
      sd = sqrt(max(sample_variance(x), 0.0_dp))
      if (sd <= tiny(1.0_dp)) then
         value = 0.0_dp
      else
         value = sum(((x - m) / sd)**4) / real(size(x), dp)
      end if
   end function sample_kurtosis

   pure real(dp) function covariance(x, y) result(value)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: mx, my
      integer :: n
      n = min(size(x), size(y))
      if (n < 2) then
         value = 0.0_dp
         return
      end if
      mx = sum(x(:n)) / real(n, dp)
      my = sum(y(:n)) / real(n, dp)
      value = sum((x(:n) - mx) * (y(:n) - my)) / real(n - 1, dp)
   end function covariance

   pure logical function all_finite(x) result(ok)
      real(dp), intent(in) :: x(:)
      integer :: i
      ok = .true.
      do i = 1, size(x)
         if (.not. ieee_is_finite(x(i))) then
            ok = .false.
            return
         end if
      end do
   end function all_finite

end module garchsk_stats
