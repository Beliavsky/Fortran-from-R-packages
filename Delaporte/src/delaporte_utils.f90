! Copyright (c) 2016, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! Derived from src/utils.f90 in the Delaporte R package.

module delaporte_utils
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use delaporte_kinds, only : dp
    implicit none
    private

    real(dp), parameter, public :: zero = 0.0_dp
    real(dp), parameter, public :: half = 0.5_dp
    real(dp), parameter, public :: one = 1.0_dp
    real(dp), parameter, public :: three_halves = 1.5_dp
    real(dp), parameter, public :: two = 2.0_dp
    real(dp), parameter, public :: three = 3.0_dp

    public :: log1p_delap, clamp01, recycle_index, nan_dp

contains

    pure elemental function log1p_delap(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        if (abs(x) <= 1.0e-4_dp) then
            y = (-half * x + one) * x
        else
            y = log(one + x)
        end if
    end function log1p_delap

    pure elemental function clamp01(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = max(zero, min(one, x))
    end function clamp01

    pure elemental function recycle_index(i, n) result(j)
        integer, intent(in) :: i, n
        integer :: j

        j = modulo(i - 1, n) + 1
    end function recycle_index

    pure function nan_dp() result(x)
        real(dp) :: x

        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_dp

end module delaporte_utils
