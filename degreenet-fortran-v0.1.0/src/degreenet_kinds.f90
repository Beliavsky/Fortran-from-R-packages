! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: huge_neg = -1.0e300_dp
end module degreenet_kinds
