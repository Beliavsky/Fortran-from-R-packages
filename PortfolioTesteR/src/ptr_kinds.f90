! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_kinds
  use, intrinsic :: iso_fortran_env, only : int64
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: i8 = int64
end module ptr_kinds
