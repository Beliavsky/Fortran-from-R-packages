! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module multiatsm_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module multiatsm_kinds
