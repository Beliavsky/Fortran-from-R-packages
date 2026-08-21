! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_rng
    use iso_fortran_env, only : int64
    use dirmult_types, only : dp
    implicit none
    private
    public :: seed_rng, normal_random, gamma_random, multinomial_random

contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer, allocatable :: put(:)
        integer :: n, i
        integer(int64) :: z, h

        call random_seed(size=n)
        allocate(put(n))
        h = int(huge(1), int64) - 1_int64
        do i = 1, n
            z = modulo(int(seed, int64) + 104729_int64 * int(i, int64) + &
                       15485863_int64 * int(i*i, int64), h)
            put(i) = int(max(1_int64, z))
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    function normal_random(mean, sd) result(x)
        real(dp), intent(in), optional :: mean, sd
        real(dp) :: x, mu, sigma, u1, u2
        real(dp), parameter :: twopi = 6.283185307179586476925286766559_dp

        mu = 0.0_dp
        sigma = 1.0_dp
        if (present(mean)) mu = mean
        if (present(sd)) sigma = sd
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        x = mu + sigma * sqrt(-2.0_dp * log(u1)) * cos(twopi * u2)
    end function normal_random

    recursive function gamma_random(shape) result(x)
        real(dp), intent(in) :: shape
        real(dp) :: x
        real(dp) :: d, c, z, v, u

        if (shape <= 0.0_dp) then
            x = 0.0_dp
            return
        end if
        if (shape < 1.0_dp) then
            call random_number(u)
            u = max(u, tiny(1.0_dp))
            x = gamma_random(shape + 1.0_dp) * u ** (1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            z = normal_random()
            v = 1.0_dp + c * z
            if (v <= 0.0_dp) cycle
            v = v * v * v
            call random_number(u)
            if (u < 1.0_dp - 0.0331_dp * z**4) exit
            if (log(max(u, tiny(1.0_dp))) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
        end do
        x = d * v
    end function gamma_random

    subroutine multinomial_random(n, p, x)
        integer, intent(in) :: n
        real(dp), intent(in) :: p(:)
        integer, intent(out) :: x(:)
        real(dp), allocatable :: cdf(:)
        real(dp) :: s, u
        integer :: k, i, j

        k = size(p)
        x = 0
        if (size(x) /= k .or. n <= 0 .or. k == 0) return
        s = sum(p)
        if (s <= 0.0_dp) return
        allocate(cdf(k))
        cdf(1) = p(1) / s
        do j = 2, k
            cdf(j) = cdf(j-1) + p(j) / s
        end do
        cdf(k) = 1.0_dp
        do i = 1, n
            call random_number(u)
            do j = 1, k
                if (u <= cdf(j)) then
                    x(j) = x(j) + 1
                    exit
                end if
            end do
        end do
    end subroutine multinomial_random

end module dirmult_rng
