! Modern Fortran translation of R package skewunit.
! SPDX-License-Identifier: GPL-2.0-or-later
module skewunit_kinds
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter, public :: sqrt2pi = sqrt(2.0_dp*pi)
   real(dp), parameter, public :: eps_dp = epsilon(1.0_dp)

   public :: nan_dp, pos_inf_dp, neg_inf_dp

contains

   pure real(dp) function nan_dp() result(x)
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure real(dp) function pos_inf_dp() result(x)
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function pos_inf_dp

   pure real(dp) function neg_inf_dp() result(x)
      x = ieee_value(0.0_dp, ieee_negative_inf)
   end function neg_inf_dp

end module skewunit_kinds
