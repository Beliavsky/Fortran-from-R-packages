! SPDX-License-Identifier: GPL-2.0-or-later
module mc2d_kinds
  use iso_fortran_env, only : real64, int64
  implicit none
  private
  public :: dp, i64, pi, nan_dp
  integer, parameter :: dp = real64
  integer, parameter :: i64 = int64
  real(dp), parameter :: pi = acos(-1.0_dp)
contains
  real(dp) function nan_dp() result(x)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    x = ieee_value(0.0_dp, ieee_quiet_nan)
  end function nan_dp
end module mc2d_kinds
