! Upstream OOR license declaration: LGPL (version unspecified).
module oor_test_functions
   use oor_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: guirland, sin1, difficult, difficult2, double_sine
contains
   elemental function guirland(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = 4.0_dp * x * (1.0_dp - x) * &
          (0.75_dp + 0.25_dp * (1.0_dp - sqrt(abs(sin(60.0_dp * x)))))
   end function guirland

   elemental function sin1(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      y = 0.5_dp * sin(13.0_dp * x) * sin(27.0_dp * x) + 0.5_dp
   end function sin1

   elemental function difficult(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) <= tiny(1.0_dp)) then
         y = ieee_value(y, ieee_quiet_nan)
      else
         y = 1.0_dp - sqrt(x) + (-x*x + sqrt(x)) * &
             (sin(1.0_dp / (x*x*x)) + 1.0_dp) / 2.0_dp
      end if
   end function difficult

   elemental function difficult2(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y, z, tmp
      z = abs(x - 0.5_dp)
      if (z <= tiny(1.0_dp)) then
         y = ieee_value(y, ieee_quiet_nan)
         return
      end if
      tmp = log(z) / log(2.0_dp)
      if (modulo(tmp, 1.0_dp) <= 0.5_dp) then
         y = sqrt(z) - (x - 0.5_dp)**2 - sqrt(z)
      else
         y = -sqrt(z)
      end if
   end function difficult2

   elemental function double_sine(x, rho1, rho2, tmax) result(y)
      real(dp), intent(in) :: x
      real(dp), intent(in), optional :: rho1, rho2, tmax
      real(dp) :: y, r1, r2, tm, u, envelope_width, phase
      real(dp), parameter :: pi = acos(-1.0_dp)

      r1 = 0.3_dp
      r2 = 0.8_dp
      tm = 0.5_dp
      if (present(rho1)) r1 = rho1
      if (present(rho2)) r2 = rho2
      if (present(tmax)) tm = tmax

      u = 2.0_dp * abs(x - tm)
      if (u <= tiny(1.0_dp)) then
         y = 0.0_dp
         return
      end if
      envelope_width = u**(-log(r2)/log(2.0_dp)) - u**(-log(r1)/log(2.0_dp))
      phase = (log(u)/log(2.0_dp)) / 2.0_dp
      y = 0.5_dp * (sin(2.0_dp * phase * pi) + 1.0_dp) * envelope_width - &
          u**(-log(r2)/log(2.0_dp))
   end function double_sine
end module oor_test_functions
