! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from tsDyn, copyright its original authors.
! This file may be redistributed and/or modified under GPL version 2 or later.
module tsdyn_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
  public :: n_deterministic
  integer, parameter, public :: include_none = 0
  integer, parameter, public :: include_const = 1
  integer, parameter, public :: include_trend = 2
  integer, parameter, public :: include_both = 3
contains
  pure integer function n_deterministic(include) result(n)
    integer, intent(in) :: include
    select case (include)
    case (include_none)
      n = 0
    case (include_const, include_trend)
      n = 1
    case (include_both)
      n = 2
    case default
      n = -1
    end select
  end function n_deterministic
end module tsdyn_kinds
