! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_kinds
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
  real(dp), parameter, public :: sqrt_two_pi = sqrt(2.0_dp*pi)
  real(dp), parameter, public :: log_two_pi = log(2.0_dp*pi)
  real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)
  real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module sn_kinds
