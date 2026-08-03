! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = selected_int_kind(18)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)
  real(dp), parameter, public :: huge_dp = huge(1.0_dp)
end module r4gpf_kinds
