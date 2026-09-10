! SPDX-License-Identifier: GPL-2.0-only
module combinat_rng
    use, intrinsic :: iso_fortran_env, only : int64
    use combinat_kinds, only : dp
    implicit none
    private

    integer(int64), parameter :: park_miller_modulus = 2147483647_int64
    integer(int64), parameter :: park_miller_multiplier = 16807_int64
    integer(int64), parameter :: park_miller_q = 127773_int64
    integer(int64), parameter :: park_miller_r = 2836_int64

    type, public :: combinat_rng_state
        integer(int64) :: state = 1234567_int64
    end type combinat_rng_state

    public :: rng_seed, rng_uniform

contains

    pure subroutine rng_seed(rng, seed)
        type(combinat_rng_state), intent(out) :: rng !! RNG state initialized for reproducible multinomial simulation.
        integer(int64), intent(in) :: seed !! User seed mapped to the Park-Miller state range.

        integer(int64) :: mapped

        mapped = modulo(seed, park_miller_modulus - 1_int64)
        if (mapped == 0_int64) then
            rng%state = 1_int64
        else
            rng%state = mapped
        end if
    end subroutine rng_seed

    real(dp) function rng_uniform(rng) result(value)
        type(combinat_rng_state), intent(inout) :: rng !! RNG state advanced by one Park-Miller uniform draw.

        integer(int64) :: quotient
        integer(int64) :: next_state

        quotient = rng%state/park_miller_q
        next_state = park_miller_multiplier*(rng%state - quotient*park_miller_q) - park_miller_r*quotient
        if (next_state <= 0_int64) next_state = next_state + park_miller_modulus
        rng%state = next_state
        value = real(next_state, dp)/real(park_miller_modulus, dp)
    end function rng_uniform

end module combinat_rng
