! SPDX-License-Identifier: GPL-2.0-only
module ks_kinds
   use iso_fortran_env, only: real64, int64
   implicit none
   private
   public :: dp, i8, pi, sqrt2, log2pi
   integer, parameter :: dp = real64
   integer, parameter :: i8 = int64
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter :: log2pi = log(2.0_dp*pi)
end module ks_kinds
