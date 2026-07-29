! FatTailsR modern Fortran translation
! Copyright (C) 2014-2026 Patrice Kiener
! Licensed under GPL-2.0-only. See COPYING.
module fattailsr_moments
   use fattailsr_kinds, only : dp
   use fattailsr_math, only : pi, sqrt3, beta_fn, binomial_real
   use fattailsr_params, only : kiener_parameters
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private

   type, public :: moment_summary
      real(dp) :: mean = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      real(dp) :: skewness = 0.0_dp
      real(dp) :: kurtosis = 0.0_dp
      real(dp) :: excess_kurtosis = 0.0_dp
   end type moment_summary

   public :: kiener_raw_moment, kiener_central_moment, kiener_moment_summary
   public :: sample_moment_summary

contains

   pure function kiener_raw_moment(n, par) result(value)
      integer, intent(in) :: n
      type(kiener_parameters), intent(in) :: par
      real(dp) :: value
      integer :: i, j
      real(dp) :: scale, term, aa, bb

      if (n < 0 .or. n >= min(par%a, par%w)) then
         value = ieee_value(value, ieee_quiet_nan)
         return
      end if
      scale = sqrt3/pi/2.0_dp*par%g*par%k
      value = 0.0_dp
      do i = 0, n
         do j = 0, i
            aa = 1.0_dp - real(j,dp)/par%a + real(i-j,dp)/par%w
            bb = 1.0_dp + real(j,dp)/par%a - real(i-j,dp)/par%w
            term = binomial_real(n,i)*par%m**(n-i)*scale**i*&
                   binomial_real(i,j)*(-1.0_dp)**j*beta_fn(aa,bb)
            value = value + term
         end do
      end do
   end function kiener_raw_moment

   pure function kiener_central_moment(n, par) result(value)
      integer, intent(in) :: n
      type(kiener_parameters), intent(in) :: par
      real(dp) :: value
      integer :: i, j
      real(dp) :: scale, nu, reduced, aa, bb

      if (n < 0 .or. n >= min(par%a, par%w)) then
         value = ieee_value(value, ieee_quiet_nan)
         return
      end if
      nu = -beta_fn(1.0_dp - 1.0_dp/par%a, 1.0_dp + 1.0_dp/par%a) + &
            beta_fn(1.0_dp - 1.0_dp/par%w, 1.0_dp + 1.0_dp/par%w)
      reduced = 0.0_dp
      do i = 0, n
         do j = 0, i
            aa = 1.0_dp - real(j,dp)/par%a + real(i-j,dp)/par%w
            bb = 1.0_dp + real(j,dp)/par%a - real(i-j,dp)/par%w
            reduced = reduced + binomial_real(n,i)*binomial_real(i,j)*&
                      (-nu)**(n-i)*(-1.0_dp)**j*beta_fn(aa,bb)
         end do
      end do
      scale = sqrt3/pi/2.0_dp*par%g*par%k
      value = reduced*scale**n
   end function kiener_central_moment

   pure function kiener_moment_summary(par) result(summary)
      type(kiener_parameters), intent(in) :: par
      type(moment_summary) :: summary
      real(dp) :: u2, u3, u4
      summary%mean = kiener_raw_moment(1, par)
      u2 = kiener_central_moment(2, par)
      u3 = kiener_central_moment(3, par)
      u4 = kiener_central_moment(4, par)
      summary%standard_deviation = sqrt(u2)
      summary%skewness = u3/u2**1.5_dp
      summary%kurtosis = u4/u2**2
      summary%excess_kurtosis = summary%kurtosis - 3.0_dp
   end function kiener_moment_summary

   pure function sample_moment_summary(x) result(summary)
      real(dp), intent(in) :: x(:)
      type(moment_summary) :: summary
      real(dp) :: mean, u2, u3, u4, n
      if (size(x) < 2) then
         summary%mean = ieee_value(mean, ieee_quiet_nan)
         summary%standard_deviation = summary%mean
         summary%skewness = summary%mean
         summary%kurtosis = summary%mean
         summary%excess_kurtosis = summary%mean
         return
      end if
      n = real(size(x),dp)
      mean = sum(x)/n
      u2 = sum((x-mean)**2)/n
      u3 = sum((x-mean)**3)/n
      u4 = sum((x-mean)**4)/n
      summary%mean = mean
      summary%standard_deviation = sqrt(sum((x-mean)**2)/real(size(x)-1,dp))
      summary%skewness = u3/summary%standard_deviation**3
      summary%kurtosis = u4/summary%standard_deviation**4
      summary%excess_kurtosis = summary%kurtosis - 3.0_dp
   end function sample_moment_summary

end module fattailsr_moments
