! SPDX-License-Identifier: MIT
module jumptest_kinds
  use, intrinsic :: iso_fortran_env, only : int64, real64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)

end module jumptest_kinds
