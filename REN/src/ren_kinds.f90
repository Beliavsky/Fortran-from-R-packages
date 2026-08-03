! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_kinds
  use iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
end module ren_kinds
