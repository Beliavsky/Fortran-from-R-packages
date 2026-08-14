! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn_special
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf
  use fnn_kinds, only : dp
  implicit none
  private
  public :: digamma_dp, log_unit_ball_volume
contains
  pure elemental real(dp) function digamma_dp(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z, inv, inv2
    if (x <= 0.0_dp) then
      y = ieee_value(0.0_dp,ieee_negative_inf)
      return
    end if
    z = x
    y = 0.0_dp
    do while (z < 8.0_dp)
      y = y - 1.0_dp/z
      z = z + 1.0_dp
    end do
    inv = 1.0_dp/z
    inv2 = inv*inv
    y = y + log(z) - 0.5_dp*inv - inv2*(1.0_dp/12.0_dp - inv2*(1.0_dp/120.0_dp &
      - inv2*(1.0_dp/252.0_dp - inv2*(1.0_dp/240.0_dp - inv2*5.0_dp/660.0_dp))))
  end function digamma_dp

  pure real(dp) function log_unit_ball_volume(p) result(v)
    integer, intent(in) :: p
    real(dp), parameter :: pi = acos(-1.0_dp)
    v = 0.5_dp*real(p,dp)*log(pi) - log_gamma(0.5_dp*real(p,dp) + 1.0_dp)
  end function log_unit_ball_volume
end module fnn_special
