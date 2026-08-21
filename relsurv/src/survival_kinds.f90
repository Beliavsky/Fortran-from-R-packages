! SPDX-License-Identifier: LGPL-2.0-or-later
module survival_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module survival_kinds
