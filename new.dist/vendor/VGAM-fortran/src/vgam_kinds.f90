! SPDX-License-Identifier: GPL-3.0-only
module vgam_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int32, int64
   implicit none
   private
   integer, parameter, public :: dp = real64
   integer, parameter, public :: i32 = int32
   integer, parameter, public :: i64 = int64
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter, public :: log2pi = log(2.0_dp*pi)
end module vgam_kinds
