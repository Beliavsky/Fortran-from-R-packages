! SPDX-License-Identifier: GPL-3.0-only
module esback_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
  real(dp), parameter, public :: huge_penalty = 1.0e100_dp
end module esback_kinds
