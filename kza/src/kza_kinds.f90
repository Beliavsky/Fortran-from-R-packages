! SPDX-License-Identifier: GPL-3.0-only
module kza_kinds
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: pi = acos(-1.0_dp)

end module kza_kinds
