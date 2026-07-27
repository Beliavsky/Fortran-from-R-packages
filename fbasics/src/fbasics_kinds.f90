! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt2 = sqrt(2.0_dp)
  real(dp), parameter, public :: tiny_dp = tiny(1.0_dp)
  public :: clamp
contains
  pure elemental real(dp) function clamp(x, lo, hi) result(y)
    real(dp), intent(in) :: x, lo, hi
    y = min(max(x, lo), hi)
  end function clamp
end module fbasics_kinds
