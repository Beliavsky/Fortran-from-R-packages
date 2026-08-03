! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = selected_int_kind(18)
  integer, parameter, public :: status_ok = 0
  integer, parameter, public :: status_invalid = 1
  integer, parameter, public :: status_singular = 2
  integer, parameter, public :: status_nonconverged = 3
  integer, parameter, public :: status_numerical = 4
end module intrinsicfrp_kinds
