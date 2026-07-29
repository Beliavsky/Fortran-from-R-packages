! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of cvar 0.6 by Georgi N. Boshnakov.
module cvar_rng
    use, intrinsic :: iso_fortran_env, only : int64
    use cvar_kinds, only : dp
    use cvar_special, only : ged_scale
    implicit none
    private

    type, public :: rng_state
        private
        integer(int64) :: state = 1_int64
        logical :: has_spare = .false.
        real(dp) :: spare = 0.0_dp
    end type rng_state

    public :: rng_seed, rng_uniform, rng_normal, rng_gamma
    public :: rng_student_t, rng_std_student_t, rng_ged

    integer(int64), parameter :: multiplier = 16807_int64
    integer(int64), parameter :: modulus = 2147483647_int64
    integer(int64), parameter :: quotient = 127773_int64
    integer(int64), parameter :: remainder = 2836_int64
    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    subroutine rng_seed(rng, seed)
        type(rng_state), intent(inout) :: rng
        integer(int64), intent(in) :: seed
        integer(int64) :: s

        s = modulo(seed, modulus - 1_int64)
        if (s < 0_int64) s = s + modulus - 1_int64
        rng%state = s + 1_int64
        rng%has_spare = .false.
        rng%spare = 0.0_dp
    end subroutine rng_seed

    function rng_uniform(rng) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp) :: value
        integer(int64) :: hi, lo, test

        hi = rng%state / quotient
        lo = modulo(rng%state, quotient)
        test = multiplier * lo - remainder * hi
        if (test > 0_int64) then
            rng%state = test
        else
            rng%state = test + modulus
        end if
        value = real(rng%state, dp) / real(modulus, dp)
        value = min(1.0_dp - epsilon(1.0_dp), max(tiny(1.0_dp), value))
    end function rng_uniform

    function rng_normal(rng) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp) :: value
        real(dp) :: radius, theta

        if (rng%has_spare) then
            value = rng%spare
            rng%has_spare = .false.
            return
        end if

        radius = sqrt(-2.0_dp * log(rng_uniform(rng)))
        theta = 2.0_dp * pi * rng_uniform(rng)
        value = radius * cos(theta)
        rng%spare = radius * sin(theta)
        rng%has_spare = .true.
    end function rng_normal

    recursive function rng_gamma(rng, shape) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: shape
        real(dp) :: value
        real(dp) :: d, c, x, v, u

        if (shape <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (shape < 1.0_dp) then
            value = rng_gamma(rng, shape + 1.0_dp) * rng_uniform(rng)**(1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            x = rng_normal(rng)
            v = 1.0_dp + c * x
            if (v <= 0.0_dp) cycle
            v = v * v * v
            u = rng_uniform(rng)
            if (u < 1.0_dp - 0.0331_dp * x**4) exit
            if (log(u) < 0.5_dp * x * x + d * (1.0_dp - v + log(v))) exit
        end do
        value = d * v
    end function rng_gamma

    function rng_student_t(rng, nu) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: nu
        real(dp) :: value, chi_square

        chi_square = 2.0_dp * rng_gamma(rng, 0.5_dp * nu)
        value = rng_normal(rng) / sqrt(chi_square / nu)
    end function rng_student_t

    function rng_std_student_t(rng, nu) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: nu
        real(dp) :: value

        value = sqrt((nu - 2.0_dp) / nu) * rng_student_t(rng, nu)
    end function rng_std_student_t

    function rng_ged(rng, nu) result(value)
        type(rng_state), intent(inout) :: rng
        real(dp), intent(in) :: nu
        real(dp) :: value, magnitude

        magnitude = (2.0_dp * rng_gamma(rng, 1.0_dp / nu))**(1.0_dp / nu)
        value = ged_scale(nu) * magnitude
        if (rng_uniform(rng) < 0.5_dp) value = -value
    end function rng_ged

end module cvar_rng
