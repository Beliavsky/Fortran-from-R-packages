! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from MCMCglmm 2.36 computational code; see NOTICE.md and upstream/.
module mcmcglmm_rng
    use iso_fortran_env, only : int64
    use r_kinds, only : dp
    implicit none
    private

    integer(int64), parameter :: pm_modulus = 2147483647_int64
    integer(int64), parameter :: pm_multiplier = 16807_int64
    integer(int64), parameter :: pm_q = 127773_int64
    integer(int64), parameter :: pm_r = 2836_int64
    real(dp), parameter :: two_pi = 6.2831853071795864769252867665590058_dp

    type, public :: rng_state
        integer(int64) :: seed = 1_int64
        logical :: has_spare = .false.
        real(dp) :: spare = 0.0_dp
    end type rng_state

    public :: rng_seed
    public :: rng_uniform
    public :: rng_normal
    public :: rng_exponential
    public :: rng_gamma
    public :: rng_chisq
    public :: rng_poisson

contains

    pure subroutine rng_seed(state, seed)
        type(rng_state), intent(out) :: state !! Generator state to initialize.
        integer(int64), intent(in) :: seed !! Integer seed; zero maps to one and negatives fold into the valid range.
        integer(int64) :: s

        s = modulo(abs(seed), pm_modulus)
        if (s == 0_int64) s = 1_int64
        state%seed = s
        state%has_spare = .false.
        state%spare = 0.0_dp
    end subroutine rng_seed

    pure subroutine rng_uniform(state, value)
        type(rng_state), intent(inout) :: state !! Generator state advanced by one Park-Miller draw.
        real(dp), intent(out) :: value !! Uniform variate strictly between zero and one.
        integer(int64) :: hi
        integer(int64) :: lo
        integer(int64) :: test

        hi = state%seed / pm_q
        lo = state%seed - hi * pm_q
        test = pm_multiplier * lo - pm_r * hi
        if (test > 0_int64) then
            state%seed = test
        else
            state%seed = test + pm_modulus
        end if
        value = real(state%seed, dp) / real(pm_modulus, dp)
    end subroutine rng_uniform

    pure subroutine rng_normal(state, value)
        type(rng_state), intent(inout) :: state !! Generator state, including the cached Box-Muller companion variate.
        real(dp), intent(out) :: value !! Standard-normal variate with mean zero and variance one.
        real(dp) :: u1
        real(dp) :: u2
        real(dp) :: radius

        if (state%has_spare) then
            value = state%spare
            state%has_spare = .false.
            return
        end if

        call rng_uniform(state, u1)
        call rng_uniform(state, u2)
        radius = sqrt(-2.0_dp * log(u1))
        value = radius * cos(two_pi * u2)
        state%spare = radius * sin(two_pi * u2)
        state%has_spare = .true.
    end subroutine rng_normal

    pure subroutine rng_exponential(state, mean_value, value)
        type(rng_state), intent(inout) :: state !! Generator state advanced by one uniform draw.
        real(dp), intent(in) :: mean_value !! Positive exponential mean; values not greater than zero return zero.
        real(dp), intent(out) :: value !! Exponential variate with the requested mean.
        real(dp) :: u

        if (mean_value <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        call rng_uniform(state, u)
        value = -mean_value * log(u)
    end subroutine rng_exponential

    pure recursive subroutine rng_gamma(state, shape, scale, value)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Marsaglia-Tsang gamma sampler.
        real(dp), intent(in) :: shape !! Positive gamma shape parameter.
        real(dp), intent(in) :: scale !! Positive gamma scale parameter.
        real(dp), intent(out) :: value !! Gamma variate; zero is returned for invalid nonpositive parameters.
        real(dp) :: d
        real(dp) :: c
        real(dp) :: x
        real(dp) :: v
        real(dp) :: u
        real(dp) :: g

        if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        if (shape < 1.0_dp) then
            call rng_gamma(state, shape + 1.0_dp, 1.0_dp, g)
            call rng_uniform(state, u)
            value = scale * g * u ** (1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            call rng_normal(state, x)
            v = 1.0_dp + c * x
            if (v <= 0.0_dp) cycle
            v = v * v * v
            call rng_uniform(state, u)
            if (u < 1.0_dp - 0.0331_dp * x ** 4) exit
            if (log(u) < 0.5_dp * x * x + d * (1.0_dp - v + log(v))) exit
        end do
        value = scale * d * v
    end subroutine rng_gamma

    pure subroutine rng_chisq(state, degrees_freedom, value)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the gamma representation of chi-square sampling.
        real(dp), intent(in) :: degrees_freedom !! Positive chi-square degrees of freedom.
        real(dp), intent(out) :: value !! Chi-square variate; zero is returned for invalid degrees of freedom.

        call rng_gamma(state, 0.5_dp * degrees_freedom, 2.0_dp, value)
    end subroutine rng_chisq

    pure subroutine rng_poisson(state, lambda, value)
        type(rng_state), intent(inout) :: state !! Generator state consumed by the Poisson sampler.
        real(dp), intent(in) :: lambda !! Nonnegative Poisson mean.
        integer, intent(out) :: value !! Nonnegative Poisson count; zero is returned when lambda is nonpositive.
        real(dp) :: threshold
        real(dp) :: product
        real(dp) :: u
        real(dp) :: z
        integer :: k

        if (lambda <= 0.0_dp) then
            value = 0
            return
        end if

        if (lambda < 30.0_dp) then
            threshold = exp(-lambda)
            product = 1.0_dp
            k = 0
            do
                k = k + 1
                call rng_uniform(state, u)
                product = product * u
                if (product <= threshold) exit
            end do
            value = k - 1
            return
        end if

        do
            call rng_normal(state, z)
            k = nint(lambda + sqrt(lambda) * z)
            if (k >= 0) exit
        end do
        value = k
    end subroutine rng_poisson

end module mcmcglmm_rng
