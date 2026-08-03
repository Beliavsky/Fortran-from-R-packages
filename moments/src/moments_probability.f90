! SPDX-License-Identifier: GPL-2.0-or-later
module moments_probability
   use moments_kinds, only : dp
   implicit none
   private
   public :: normal_cdf, normal_survival, chi_square_2_survival

contains

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function normal_survival(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(x / sqrt(2.0_dp))
   end function normal_survival

   elemental real(dp) function chi_square_2_survival(x) result(p)
      real(dp), intent(in) :: x
      if (x <= 0.0_dp) then
         p = 1.0_dp
      else
         p = exp(-0.5_dp * x)
      end if
   end function chi_square_2_survival

end module moments_probability
