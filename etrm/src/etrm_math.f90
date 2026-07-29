! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_math
   use etrm_kinds, only : dp
   implicit none
   private

   public :: normal_cdf, signum, trade_round

contains

   elemental real(dp) function normal_cdf(x) result(p)
      real(dp), intent(in) :: x
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   elemental real(dp) function signum(x) result(s)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         s = 1.0_dp
      else if (x < 0.0_dp) then
         s = -1.0_dp
      else
         s = 0.0_dp
      end if
   end function signum

   elemental real(dp) function trade_round(x, integer_trades) result(y)
      real(dp), intent(in) :: x
      logical, intent(in) :: integer_trades
      real(dp) :: ax, base, frac, rounded
      integer :: ibase

      if (.not. integer_trades) then
         y = x
         return
      end if

      ax = abs(x)
      base = floor(ax)
      frac = ax - base
      ibase = int(base)

      if (frac < 0.5_dp - 8.0_dp * epsilon(ax)) then
         rounded = base
      else if (frac > 0.5_dp + 8.0_dp * epsilon(ax)) then
         rounded = base + 1.0_dp
      else if (mod(ibase, 2) == 0) then
         rounded = base
      else
         rounded = base + 1.0_dp
      end if

      y = sign(rounded, x)
   end function trade_round

end module etrm_math
