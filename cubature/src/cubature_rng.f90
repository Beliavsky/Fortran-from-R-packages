! SPDX-License-Identifier: GPL-3.0-or-later
module cubature_rng
    use cubature_kinds, only : dp, i8
    implicit none
    private
    public :: rng_state, rng_seed, rng_uniform, rng_normal

    type :: rng_state
        integer(i8) :: state = 1_i8
        logical :: have_normal = .false.
        real(dp) :: normal_cache = 0.0_dp
    end type rng_state

contains

    subroutine rng_seed(rng, seed)
        type(rng_state), intent(inout) :: rng
        integer, intent(in) :: seed
        integer(i8) :: s
        s = modulo(int(seed, i8), 2147483646_i8)
        if (s <= 0_i8) s = s + 2147483646_i8
        rng%state = s
        rng%have_normal = .false.
    end subroutine rng_seed

    real(dp) function rng_uniform(rng)
        type(rng_state), intent(inout) :: rng
        integer(i8) :: hi, lo, test
        ! Park-Miller minimal standard, Schrage form: no integer overflow.
        hi = rng%state / 127773_i8
        lo = modulo(rng%state, 127773_i8)
        test = 16807_i8 * lo - 2836_i8 * hi
        if (test <= 0_i8) test = test + 2147483647_i8
        rng%state = test
        rng_uniform = real(test, dp) / 2147483647.0_dp
    end function rng_uniform

    real(dp) function rng_normal(rng)
        type(rng_state), intent(inout) :: rng
        real(dp) :: u1, u2, r, theta
        real(dp), parameter :: twopi = 6.28318530717958647692528676655900577_dp
        if (rng%have_normal) then
            rng_normal = rng%normal_cache
            rng%have_normal = .false.
            return
        end if
        u1 = max(rng_uniform(rng), tiny(1.0_dp))
        u2 = rng_uniform(rng)
        r = sqrt(-2.0_dp * log(u1))
        theta = twopi * u2
        rng_normal = r * cos(theta)
        rng%normal_cache = r * sin(theta)
        rng%have_normal = .true.
    end function rng_normal

end module cubature_rng
