! SPDX-License-Identifier: GPL-3.0-only
module nmof_kinds
   use, intrinsic :: iso_fortran_env, only: int64, real64
   implicit none
   private
   public :: dp, i8, pi
   integer, parameter :: dp = real64
   integer, parameter :: i8 = int64
   real(dp), parameter :: pi = acos(-1.0_dp)
end module nmof_kinds
