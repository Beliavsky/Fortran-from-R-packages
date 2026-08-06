! SPDX-License-Identifier: GPL-3.0-only
module spantest_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi_dp = acos(-1.0_dp)
end module spantest_kinds
