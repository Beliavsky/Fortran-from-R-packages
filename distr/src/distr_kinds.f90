! distr-fortran -- computational translation of the R package distr.
! Copyright (C) 2005-2025 distr authors.
! SPDX-License-Identifier: LGPL-3.0-only
module distr_kinds
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private
   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter, public :: sqrt2pi = sqrt(2.0_dp*pi)
   real(dp), parameter, public :: eps_dp = epsilon(1.0_dp)
   public :: nan_dp, inf_dp
contains
   pure real(dp) function nan_dp() result(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure real(dp) function inf_dp(sign) result(x)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf
      integer, intent(in), optional :: sign
      if (present(sign)) then
         if (sign < 0) then
            x = ieee_value(0.0_dp, ieee_negative_inf)
         else
            x = ieee_value(0.0_dp, ieee_positive_inf)
         end if
      else
         x = ieee_value(0.0_dp, ieee_positive_inf)
      end if
   end function inf_dp
end module distr_kinds
