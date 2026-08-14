! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_kinds
  use iso_fortran_env, only : real64, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
end module gpareto_kinds
