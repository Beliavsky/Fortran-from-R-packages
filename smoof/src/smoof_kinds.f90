! SPDX-License-Identifier: BSD-2-Clause
module smoof_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module smoof_kinds
