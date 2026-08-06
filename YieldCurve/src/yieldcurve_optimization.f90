! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve_optimization
   use yieldcurve_kinds, only : dp
   use yieldcurve_factors, only : factor_beta2, beta2_spot
   implicit none
   private

   public :: maximize_factor_beta2, maximize_beta2_spot

contains

   pure function maximize_factor_beta2(target_maturity, lower, upper) result(maximizer)
      real(dp), intent(in) :: target_maturity, lower, upper
      real(dp) :: maximizer

      maximizer = golden_factor(target_maturity, lower, upper)
   end function maximize_factor_beta2

   pure function maximize_beta2_spot(target_maturity, lower, upper) result(maximizer)
      real(dp), intent(in) :: target_maturity, lower, upper
      real(dp) :: maximizer

      maximizer = golden_spot(target_maturity, lower, upper)
   end function maximize_beta2_spot

   pure function golden_factor(target, lower, upper) result(xbest)
      real(dp), intent(in) :: target, lower, upper
      real(dp) :: xbest
      real(dp) :: a, b, c, d, fc, fd, tolerance
      real(dp), parameter :: invphi = 0.6180339887498948482_dp
      integer :: iter

      a = lower
      b = upper
      tolerance = epsilon(1.0_dp)**0.25_dp
      c = b - invphi * (b - a)
      d = a + invphi * (b - a)
      fc = factor_beta2(c, target)
      fd = factor_beta2(d, target)

      do iter = 1, 200
         if (abs(b - a) <= tolerance * (abs(c) + abs(d) + 1.0_dp)) exit
         if (fc > fd) then
            b = d
            d = c
            fd = fc
            c = b - invphi * (b - a)
            fc = factor_beta2(c, target)
         else
            a = c
            c = d
            fc = fd
            d = a + invphi * (b - a)
            fd = factor_beta2(d, target)
         end if
      end do
      xbest = 0.5_dp * (a + b)
   end function golden_factor

   pure function golden_spot(target, lower, upper) result(xbest)
      real(dp), intent(in) :: target, lower, upper
      real(dp) :: xbest
      real(dp) :: a, b, c, d, fc, fd, tolerance
      real(dp), parameter :: invphi = 0.6180339887498948482_dp
      integer :: iter

      a = lower
      b = upper
      tolerance = epsilon(1.0_dp)**0.25_dp
      c = b - invphi * (b - a)
      d = a + invphi * (b - a)
      fc = beta2_spot(target, c)
      fd = beta2_spot(target, d)

      do iter = 1, 200
         if (abs(b - a) <= tolerance * (abs(c) + abs(d) + 1.0_dp)) exit
         if (fc > fd) then
            b = d
            d = c
            fd = fc
            c = b - invphi * (b - a)
            fc = beta2_spot(target, c)
         else
            a = c
            c = d
            fc = fd
            d = a + invphi * (b - a)
            fd = beta2_spot(target, d)
         end if
      end do
      xbest = 0.5_dp * (a + b)
   end function golden_spot

end module yieldcurve_optimization
