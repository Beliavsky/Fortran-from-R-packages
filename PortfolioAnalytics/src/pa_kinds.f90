! SPDX-License-Identifier: GPL-3.0-only
module pa_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = selected_int_kind(18)
  real(dp), parameter, public :: pa_huge = huge(1.0_dp) / 1000.0_dp
  real(dp), parameter, public :: pa_pi = acos(-1.0_dp)
end module pa_kinds
