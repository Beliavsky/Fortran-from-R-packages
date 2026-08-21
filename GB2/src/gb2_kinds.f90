! SPDX-License-Identifier: GPL-2.0-or-later
module gb2_kinds
  use, intrinsic :: iso_fortran_env, only : real64
  implicit none
  integer, parameter, public :: dp = real64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module gb2_kinds
