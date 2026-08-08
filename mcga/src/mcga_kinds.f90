! SPDX-License-Identifier: GPL-2.0-or-later
module mcga_kinds
  use, intrinsic :: iso_fortran_env, only : real64, int8, int32
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int8
  integer, parameter, public :: i32 = int32
end module mcga_kinds
