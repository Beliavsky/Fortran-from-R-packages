! SPDX-License-Identifier: GPL-3.0-only
! Native Fortran translation of the computational core of rvinecopulib.
module rvine_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
  real(dp), parameter, public :: eps_prob = 1.0e-12_dp
  public :: clamp_prob
contains
  pure elemental real(dp) function clamp_prob(x) result(y)
    real(dp), intent(in) :: x
    y = min(1.0_dp - eps_prob, max(eps_prob, x))
  end function clamp_prob
end module rvine_kinds
