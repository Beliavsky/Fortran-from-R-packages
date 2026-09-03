! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from R package mice 3.19.0 by Stef van Buuren, Karin Groothuis-Oudshoorn,
! and mice contributors; see NOTICE.md and PROVENANCE.md for attribution.
! Computational translation derived from mice 3.19.0.
module mice_rng
    use iso_fortran_env, only : int64
    use r_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp
    real(dp), parameter :: inv_two53 = 1.1102230246251565404236316680908203125e-16_dp

    type, public :: mice_rng_state
        integer(int64) :: state = 88172645463393265_int64
        logical :: have_spare = .false.
        real(dp) :: spare = 0.0_dp
    end type mice_rng_state

    public :: rng_seed
    public :: rng_uniform
    public :: rng_normal
    public :: rng_gamma
    public :: rng_chisq
    public :: rng_integer
    public :: rng_shuffle

contains

    subroutine rng_seed(rng, seed)
        type(mice_rng_state), intent(out) :: rng !! RNG state initialized from `seed`.
        integer(int64), intent(in), value :: seed !! Integer seed; zero is mapped to a fixed nonzero state.

        if (seed == 0_int64) then
            rng%state = 88172645463393265_int64
        else
            rng%state = seed
        end if
        rng%have_spare = .false.
        rng%spare = 0.0_dp
    end subroutine rng_seed

    function rng_uniform(rng) result(u)
        type(mice_rng_state), intent(inout) :: rng !! RNG state advanced by one uniform draw.
        real(dp) :: u
        integer(int64) :: x

        x = rng%state
        x = ieor(x, shiftl(x, 13))
        x = ieor(x, shiftr(x, 7))
        x = ieor(x, shiftl(x, 17))
        rng%state = x
        u = real(iand(shiftr(x, 11), int(z'001FFFFFFFFFFFFF', int64)), dp) * inv_two53
        if (u <= 0.0_dp) u = inv_two53
        if (u >= 1.0_dp) u = 1.0_dp - inv_two53
    end function rng_uniform

    function rng_normal(rng) result(z)
        type(mice_rng_state), intent(inout) :: rng !! RNG state advanced by a standard-normal draw.
        real(dp) :: z
        real(dp) :: radius, theta, u1, u2

        if (rng%have_spare) then
            z = rng%spare
            rng%have_spare = .false.
            return
        end if
        u1 = rng_uniform(rng)
        u2 = rng_uniform(rng)
        radius = sqrt(-2.0_dp * log(u1))
        theta = two_pi * u2
        z = radius * cos(theta)
        rng%spare = radius * sin(theta)
        rng%have_spare = .true.
    end function rng_normal

    recursive function rng_gamma(rng, shape) result(x)
        type(mice_rng_state), intent(inout) :: rng !! RNG state advanced by a gamma draw.
        real(dp), intent(in), value :: shape !! Gamma shape parameter; must be strictly positive.
        real(dp) :: x
        real(dp) :: c, d, u, v, z

        if (shape <= 0.0_dp) then
            x = 0.0_dp
            return
        end if
        if (shape < 1.0_dp) then
            u = rng_uniform(rng)
            x = rng_gamma(rng, shape + 1.0_dp) * u**(1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            z = rng_normal(rng)
            v = 1.0_dp + c * z
            if (v <= 0.0_dp) cycle
            v = v * v * v
            u = rng_uniform(rng)
            if (u < 1.0_dp - 0.0331_dp * z**4) exit
            if (log(u) < 0.5_dp * z * z + d * (1.0_dp - v + log(v))) exit
        end do
        x = d * v
    end function rng_gamma

    function rng_chisq(rng, df) result(x)
        type(mice_rng_state), intent(inout) :: rng !! RNG state advanced by a chi-square draw.
        real(dp), intent(in), value :: df !! Chi-square degrees of freedom; must be strictly positive.
        real(dp) :: x

        x = 2.0_dp * rng_gamma(rng, 0.5_dp * df)
    end function rng_chisq

    function rng_integer(rng, lower, upper) result(value)
        type(mice_rng_state), intent(inout) :: rng !! RNG state advanced by an integer draw.
        integer, intent(in), value :: lower !! Inclusive lower integer bound.
        integer, intent(in), value :: upper !! Inclusive upper integer bound; must be at least `lower`.
        integer :: value
        integer :: span

        span = upper - lower + 1
        if (span <= 1) then
            value = lower
        else
            value = lower + min(span - 1, int(rng_uniform(rng) * real(span, dp)))
        end if
    end function rng_integer

    subroutine rng_shuffle(rng, index)
        type(mice_rng_state), intent(inout) :: rng !! RNG state used to shuffle the index vector.
        integer, intent(inout) :: index(:) !! Integer vector shuffled in place using Fisher-Yates.
        integer :: i, j, tmp

        do i = size(index), 2, -1
            j = rng_integer(rng, 1, i)
            tmp = index(i)
            index(i) = index(j)
            index(j) = tmp
        end do
    end subroutine rng_shuffle

end module mice_rng
