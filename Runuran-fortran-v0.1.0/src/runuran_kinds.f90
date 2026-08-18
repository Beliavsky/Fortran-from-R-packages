! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Runuran 0.41 / UNU.RAN by Wolfgang Hoermann and Josef Leydold.
module runuran_kinds
  use, intrinsic :: iso_fortran_env, only : real64, int64
  implicit none
  private
  integer, parameter, public :: dp = real64
  integer, parameter, public :: i8 = int64
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
  real(dp), parameter, public :: sqrt2pi = sqrt(2.0_dp*pi)
  real(dp), parameter, public :: huge_dp = huge(1.0_dp)
  real(dp), parameter, public :: eps_dp = epsilon(1.0_dp)
  public :: clamp
contains
  pure real(dp) function clamp(x,a,b) result(y)
    real(dp), intent(in) :: x,a,b
    y = min(max(x,a),b)
  end function clamp
end module runuran_kinds
