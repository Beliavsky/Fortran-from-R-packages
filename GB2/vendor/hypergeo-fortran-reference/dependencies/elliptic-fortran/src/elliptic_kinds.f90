! SPDX-License-Identifier: GPL-2.0-only
! Modern Fortran translation of computational code from R package elliptic.
module elliptic_kinds
    use iso_fortran_env, only : real64, int64
    implicit none
    private
    integer, parameter, public :: dp = real64
    integer, parameter, public :: i8 = int64
    real(dp), parameter, public :: pi = acos(-1.0_dp)
    complex(dp), parameter, public :: ci = (0.0_dp, 1.0_dp)
    public :: close_complex, finite_complex
contains
    pure logical function close_complex(a, b, tol) result(ok)
        complex(dp), intent(in) :: a, b
        real(dp), intent(in), optional :: tol
        real(dp) :: t
        t = 8.0_dp * epsilon(1.0_dp)
        if (present(tol)) t = tol
        ok = abs(a - b) <= t * max(1.0_dp, abs(a), abs(b))
    end function close_complex

    pure logical function finite_complex(z) result(ok)
        use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
        complex(dp), intent(in) :: z
        ok = ieee_is_finite(real(z,dp)) .and. ieee_is_finite(aimag(z))
    end function finite_complex
end module elliptic_kinds
