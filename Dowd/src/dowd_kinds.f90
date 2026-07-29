! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two_pi = sqrt(2.0_dp*pi)
end module dowd_kinds
