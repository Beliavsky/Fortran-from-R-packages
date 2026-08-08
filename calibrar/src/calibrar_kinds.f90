! SPDX-License-Identifier: GPL-2.0-only
module calibrar_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi_dp = 3.141592653589793238462643383279502884_dp
end module calibrar_kinds
