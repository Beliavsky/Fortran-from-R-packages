! sna-fortran: computational translation of the R sna package.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_kinds
    use, intrinsic :: iso_fortran_env, only : real64, int64
    use, intrinsic :: ieee_arithmetic
    implicit none
    private

    integer, parameter, public :: dp = real64
    integer, parameter, public :: i8 = int64
    real(dp), parameter, public :: sna_eps = 1.0e-12_dp

    public :: sna_nan, sna_inf, is_missing, is_finite_number

contains

    pure real(dp) function sna_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function sna_nan

    pure real(dp) function sna_inf() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function sna_inf

    elemental logical function is_missing(x) result(ans)
        real(dp), intent(in) :: x
        ans = ieee_is_nan(x)
    end function is_missing

    elemental logical function is_finite_number(x) result(ans)
        real(dp), intent(in) :: x
        ans = ieee_is_finite(x)
    end function is_finite_number

end module sna_kinds
