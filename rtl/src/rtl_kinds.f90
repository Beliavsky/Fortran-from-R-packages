! SPDX-License-Identifier: MIT
! Copyright (c) 2020 RTL Authors
module rtl_kinds
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
  real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)

end module rtl_kinds
