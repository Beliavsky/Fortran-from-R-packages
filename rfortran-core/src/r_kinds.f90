! SPDX-License-Identifier: MIT
module r_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int64
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter, public :: i64 = int64
   real(dp), parameter, public :: r_pi = acos(-1.0_dp)

end module r_kinds
