! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! Portions derived from DEoptimR by Eduardo L. T. Conceicao and contributors.
! This program is free software: you may redistribute it and/or modify it
! under the terms of GNU GPL version 2, or (at your option) any later version.
module test_benchmarks
    use deoptimr_kinds, only: dp
    implicit none
    private
    public :: sphere, quadratic, equality_constraint, square1, lower_one_constraint
    public :: aluffi, becker_lago, outside_unit_constraint, unit_equality_constraint
contains
    function sphere(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum(x*x)
    end function sphere

    function quadratic(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum(x*x)
    end function quadratic

    subroutine equality_constraint(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        values(1) = x(1) + x(2) - 1.0_dp
    end subroutine equality_constraint

    function square1(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = x(1)*x(1)
    end function square1

    subroutine lower_one_constraint(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        values(1) = 1.0_dp - x(1)
    end subroutine lower_one_constraint

    function aluffi(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = 0.25_dp*x(1)**4 - 0.5_dp*x(1)**2 + 0.1_dp*x(1) + 0.5_dp*x(2)**2
    end function aluffi

    function becker_lago(x) result(value)
        real(dp), intent(in) :: x(:)
        real(dp) :: value
        value = sum((abs(x) - 5.0_dp)**2)
    end function becker_lago

    subroutine outside_unit_constraint(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        values(1) = 1.0_dp - x(1)*x(1)
    end subroutine outside_unit_constraint

    subroutine unit_equality_constraint(x, values)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: values(:)
        values(1) = x(1)*x(1) - 1.0_dp
    end subroutine unit_equality_constraint
end module test_benchmarks
