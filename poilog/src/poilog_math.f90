! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use poilog_kinds, only : dp
   implicit none
   private
   public :: pi_dp, normal_pdf, poisson_pmf, log_poisson_pmf, safe_exp, is_finite_dp

   real(dp), parameter :: pi_dp = acos(-1.0_dp)

contains

   pure elemental logical function is_finite_dp(x) result(ok)
      real(dp), intent(in) :: x
      ok = ieee_is_finite(x)
   end function is_finite_dp

   pure elemental real(dp) function safe_exp(x) result(y)
      real(dp), intent(in) :: x
      real(dp), parameter :: log_huge = log(huge(1.0_dp))
      real(dp), parameter :: log_tiny = log(tiny(1.0_dp))
      if (x >= log_huge) then
         y = huge(1.0_dp)
      else if (x <= log_tiny) then
         y = 0.0_dp
      else
         y = exp(x)
      end if
   end function safe_exp

   pure elemental real(dp) function normal_pdf(x) result(p)
      real(dp), intent(in) :: x
      p = exp(-0.5_dp*x*x) / sqrt(2.0_dp*pi_dp)
   end function normal_pdf

   pure elemental real(dp) function log_poisson_pmf(k, lambda) result(lp)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      if (k < 0 .or. lambda < 0.0_dp) then
         lp = -huge(1.0_dp)
      else if (lambda <= 0.0_dp) then
         if (k == 0) then
            lp = 0.0_dp
         else
            lp = -huge(1.0_dp)
         end if
      else
         lp = real(k,dp)*log(lambda) - lambda - log_gamma(real(k+1,dp))
      end if
   end function log_poisson_pmf

   pure elemental real(dp) function poisson_pmf(k, lambda) result(p)
      integer, intent(in) :: k
      real(dp), intent(in) :: lambda
      real(dp) :: lp
      lp = log_poisson_pmf(k, lambda)
      if (lp <= log(tiny(1.0_dp))) then
         p = 0.0_dp
      else
         p = exp(lp)
      end if
   end function poisson_pmf

end module poilog_math
