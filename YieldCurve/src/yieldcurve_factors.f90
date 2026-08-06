! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve_factors
   use yieldcurve_kinds, only : dp
   implicit none
   private

   public :: beta1_spot, beta2_spot
   public :: beta1_forward, beta2_forward
   public :: factor_beta1, factor_beta2

contains

   elemental pure function beta1_spot(maturity, tau) result(value)
      real(dp), intent(in) :: maturity, tau
      real(dp) :: value
      real(dp) :: x

      x = maturity / tau
      value = loading_one(x)
   end function beta1_spot

   elemental pure function beta2_spot(maturity, tau) result(value)
      real(dp), intent(in) :: maturity, tau
      real(dp) :: value
      real(dp) :: x

      x = maturity / tau
      value = loading_two(x)
   end function beta2_spot

   elemental pure function beta1_forward(maturity, tau) result(value)
      real(dp), intent(in) :: maturity, tau
      real(dp) :: value
      real(dp) :: x

      x = maturity / tau
      value = exp(-x)
   end function beta1_forward

   elemental pure function beta2_forward(maturity, tau) result(value)
      real(dp), intent(in) :: maturity, tau
      real(dp) :: value
      real(dp) :: x

      x = maturity / tau
      value = x * exp(-x)
   end function beta2_forward

   elemental pure function factor_beta1(lambda, maturity) result(value)
      real(dp), intent(in) :: lambda, maturity
      real(dp) :: value

      value = loading_one(lambda * maturity)
   end function factor_beta1

   elemental pure function factor_beta2(lambda, maturity) result(value)
      real(dp), intent(in) :: lambda, maturity
      real(dp) :: value

      value = loading_two(lambda * maturity)
   end function factor_beta2

   elemental pure function loading_one(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      if (abs(x) < 1.0e-5_dp) then
         value = 1.0_dp - x / 2.0_dp + x**2 / 6.0_dp - x**3 / 24.0_dp + &
            x**4 / 120.0_dp - x**5 / 720.0_dp
      else
         value = (1.0_dp - exp(-x)) / x
      end if
   end function loading_one

   elemental pure function loading_two(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value

      if (abs(x) < 1.0e-5_dp) then
         value = x / 2.0_dp - x**2 / 3.0_dp + x**3 / 8.0_dp - &
            x**4 / 30.0_dp + x**5 / 144.0_dp
      else
         value = (1.0_dp - exp(-x)) / x - exp(-x)
      end if
   end function loading_two

end module yieldcurve_factors
