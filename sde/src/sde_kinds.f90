! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
module sde_kinds
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i32 = int32
   integer, parameter, public :: i64 = int64
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
   real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)
   real(dp), parameter, public :: huge_dp = huge(1.0_dp)

end module sde_kinds
