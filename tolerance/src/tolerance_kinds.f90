! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
  real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module tolerance_kinds
