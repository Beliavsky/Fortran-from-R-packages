! SPDX-License-Identifier: GPL-3.0-only
! Derived from sharpeRratio 1.4.3 by Damien Challet.
module sharpe_rratio_statistics
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
   use ghyp_kinds, only : dp
   use sharpe_rratio_calibration, only : calibration_a, calibration_a_medium, calibration_f
   implicit none
   private

   public :: test_n, a_full, f_full, correction_b, theta_snr
   public :: sample_variance, quantile_type7

contains

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, average

      if (size(x) < 2) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      average = sum(x)/real(size(x),dp)
      value = sum((x-average)**2)/real(size(x)-1,dp)
   end function sample_variance

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: key
      integer :: i, j

      do i = 2, size(x)
         key = x(i)
         j = i-1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j+1) = x(j)
            j = j-1
         end do
         x(j+1) = key
      end do
   end subroutine sort_real

   pure function quantile_type7(x, probability) result(value)
      real(dp), intent(in) :: x(:), probability
      real(dp) :: value, h, fraction
      real(dp), allocatable :: sorted(:)
      integer :: lo, hi

      if (size(x) == 0) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      allocate(sorted(size(x)))
      sorted = x
      call sort_real(sorted)
      if (size(sorted) == 1) then
         value = sorted(1)
         return
      end if
      h = 1.0_dp+real(size(sorted)-1,dp)*min(1.0_dp,max(0.0_dp,probability))
      lo = floor(h)
      hi = ceiling(h)
      fraction = h-real(lo,dp)
      value = (1.0_dp-fraction)*sorted(lo)+fraction*sorted(hi)
   end function quantile_type7

   function test_n(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, q1, q2, variance_all
      real(dp), allocatable :: finite_x(:), low(:), middle(:), high(:)
      integer :: n, nl, nm, nh

      n = count(ieee_is_finite(x))
      if (n < 5) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      allocate(finite_x(n))
      finite_x = pack(x,ieee_is_finite(x))
      q1 = quantile_type7(finite_x,0.2_dp)
      q2 = quantile_type7(finite_x,0.8_dp)
      nl = count(finite_x <= q1)
      nm = count(finite_x > q1 .and. finite_x < q2)
      nh = count(finite_x >= q2)
      if (min(nl,nm,nh) < 2) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      allocate(low(nl),middle(nm),high(nh))
      low = pack(finite_x,finite_x <= q1)
      middle = pack(finite_x,finite_x > q1 .and. finite_x < q2)
      high = pack(finite_x,finite_x >= q2)
      variance_all = sample_variance(finite_x)
      if (.not. ieee_is_finite(variance_all) .or. variance_all <= 0.0_dp) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      value = (sample_variance(low)+sample_variance(high)- &
         2.0_dp*sample_variance(middle))*sqrt(real(n,dp))/(1.8_dp*variance_all)
   end function test_n

   pure function a_full(r0) result(value)
      real(dp), intent(in) :: r0
      real(dp) :: value

      if (r0 < 0.7_dp) then
         value = calibration_a(r0)
      else if (r0 <= 0.997_dp) then
         value = calibration_a_medium(r0)
      else
         value = exp(2.877_dp)/sqrt(252.0_dp)/(1.0_dp-r0)**0.163_dp
      end if
   end function a_full

   pure function f_full(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      if (x < 3.0_dp) then
         value = calibration_f(x)
      else
         value = calibration_f(3.0_dp)
      end if
   end function f_full

   pure function correction_b(r0, n) result(value)
      real(dp), intent(in) :: r0
      integer, intent(in) :: n
      real(dp) :: value, scaled_r0

      scaled_r0 = r0*real(n,dp)**0.42_dp
      value = exp(f_full(scaled_r0))*r0**1.6_dp
   end function correction_b

   pure function theta_snr(r0, n, nu, nu_fixed) result(value)
      real(dp), intent(in) :: r0, nu
      integer, intent(in) :: n
      logical, intent(in) :: nu_fixed
      real(dp) :: value, estimate, correction

      if (n < 1 .or. nu <= 0.0_dp .or. abs(r0) >= 1.0_dp) then
         value = ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      estimate = a_full(abs(r0))-correction_b(abs(r0),n)*nu**(-1.5_dp)
      if (estimate <= 0.0_dp) estimate = 0.0_dp
      estimate = estimate*sign(1.0_dp,r0)

      if (nu_fixed) then
         correction = 1.0_dp-(exp(3.7726_dp-6.2661_dp*log(nu))- &
            exp(-0.36538_dp*nu-1.58686_dp))
         correction = correction*(1.0_dp-0.009248_dp)
      else
         correction = 1.0_dp-(exp(4.4635_dp-6.9240_dp*log(nu))- &
            exp(-0.18441_dp*nu-3.18909_dp))
         correction = correction*(1.0_dp+exp(-3.2_dp*log(nu)+0.9_dp)-0.009248_dp)
      end if
      value = estimate*correction
   end function theta_snr

end module sharpe_rratio_statistics
