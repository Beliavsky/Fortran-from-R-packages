! SPDX-License-Identifier: GPL-3.0-only
module svdnf_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module svdnf_kinds
