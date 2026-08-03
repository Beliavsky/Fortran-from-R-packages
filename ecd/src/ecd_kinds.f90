! SPDX-License-Identifier: Artistic-2.0
module ecd_kinds
  use, intrinsic :: iso_fortran_env, only : real64, int64
  implicit none
  private
  public :: dp, i8, pi, sqrt_pi, huge_dp, ecd_ok, ecd_invalid, ecd_no_convergence
  integer, parameter :: dp = real64
  integer, parameter :: i8 = int64
  real(dp), parameter :: pi = acos(-1.0_dp)
  real(dp), parameter :: sqrt_pi = sqrt(pi)
  real(dp), parameter :: huge_dp = huge(1.0_dp)
  integer, parameter :: ecd_ok = 0
  integer, parameter :: ecd_invalid = 1
  integer, parameter :: ecd_no_convergence = 2
end module ecd_kinds
