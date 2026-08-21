! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_kinds
  use, intrinsic :: iso_fortran_env, only : real64, int32, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i32 = int32
  integer, parameter, public :: i64 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module survey_kinds
