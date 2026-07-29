! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
  real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module actuar_kinds
