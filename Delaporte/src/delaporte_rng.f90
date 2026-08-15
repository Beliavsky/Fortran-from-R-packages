! Copyright (c) 2013, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! The upstream exact random generator uses R's uniform RNG followed by the
! Delaporte quantile. This standalone version uses Fortran RANDOM_NUMBER.
! The mixture generator implements the upstream R exact=FALSE algorithm.

module delaporte_rng
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, &
        ieee_value
    use delaporte_kinds, only : dp
    use delaporte_distribution, only : qdelap
    use delaporte_utils, only : zero, one, recycle_index
    implicit none
    private

    real(dp), parameter :: pi = acos(-one)

    public :: rdelap, rdelap_vec, seed_delaporte
    public :: gamma_random, poisson_random

contains

    subroutine seed_delaporte(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        integer(kind=8) :: x

        call random_seed(size=n)
        allocate(put(n))
        x = int(seed, kind=8)
        if (x == 0_8) x = 104729_8
        do i = 1, n
            x = modulo(1664525_8 * x + 1013904223_8, 2147483647_8)
            put(i) = int(max(1_8, x))
        end do
        call random_seed(put=put)
    end subroutine seed_delaporte

    subroutine rdelap(n, alpha, beta, lambda, values, exact)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha, beta, lambda
        real(dp), intent(out) :: values(n)
        logical, intent(in), optional :: exact
        integer :: i
        real(dp) :: u, g
        logical :: use_exact

        use_exact = .true.
        if (present(exact)) use_exact = exact

        if (alpha <= zero .or. beta <= zero .or. lambda <= zero .or. &
            ieee_is_nan(alpha + beta + lambda)) then
            values = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (use_exact) then
            do i = 1, n
                call random_number(u)
                values(i) = qdelap(u, alpha, beta, lambda)
            end do
        else
            do i = 1, n
                g = gamma_random(alpha, beta)
                values(i) = real(poisson_random(g + lambda), dp)
            end do
        end if
    end subroutine rdelap

    subroutine rdelap_vec(n, alpha, beta, lambda, values, exact)
        integer, intent(in) :: n
        real(dp), intent(in) :: alpha(:), beta(:), lambda(:)
        real(dp), intent(out) :: values(n)
        logical, intent(in), optional :: exact
        integer :: i, ia, ib, il
        real(dp) :: a, b, l, u, g
        logical :: use_exact

        use_exact = .true.
        if (present(exact)) use_exact = exact
        do i = 1, n
            ia = recycle_index(i, size(alpha))
            ib = recycle_index(i, size(beta))
            il = recycle_index(i, size(lambda))
            a = alpha(ia)
            b = beta(ib)
            l = lambda(il)
            if (a <= zero .or. b <= zero .or. l <= zero .or. &
                ieee_is_nan(a + b + l)) then
                values(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            else if (use_exact) then
                call random_number(u)
                values(i) = qdelap(u, a, b, l)
            else
                g = gamma_random(a, b)
                values(i) = real(poisson_random(g + l), dp)
            end if
        end do
    end subroutine rdelap_vec

    function gamma_random(shape, scale) result(x)
        real(dp), intent(in) :: shape, scale
        real(dp) :: x
        real(dp) :: d, c, u, v, z

        if (shape <= zero .or. scale <= zero) then
            x = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (shape < one) then
            call random_number(u)
            x = gamma_random_mt(shape + one) * u ** (one / shape) * scale
            return
        end if

        d = shape - one / 3.0_dp
        c = one / sqrt(9.0_dp * d)
        do
            z = normal_random()
            v = one + c * z
            if (v <= zero) cycle
            v = v * v * v
            call random_number(u)
            if (u < one - 0.0331_dp * z ** 4) exit
            if (log(u) < 0.5_dp * z * z + d * (one - v + log(v))) exit
        end do
        x = scale * d * v
    end function gamma_random

    function gamma_random_mt(shape) result(x)
        real(dp), intent(in) :: shape
        real(dp) :: x
        real(dp) :: d, c, u, v, z

        d = shape - one / 3.0_dp
        c = one / sqrt(9.0_dp * d)
        do
            z = normal_random()
            v = one + c * z
            if (v <= zero) cycle
            v = v * v * v
            call random_number(u)
            if (u < one - 0.0331_dp * z ** 4) exit
            if (log(u) < 0.5_dp * z * z + d * (one - v + log(v))) exit
        end do
        x = d * v
    end function gamma_random_mt

    function normal_random() result(z)
        real(dp) :: z
        real(dp) :: u1, u2

        do
            call random_number(u1)
            if (u1 > zero) exit
        end do
        call random_number(u2)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
    end function normal_random

    function poisson_random(mu) result(k)
        real(dp), intent(in) :: mu
        integer :: k
        real(dp) :: a, b, inv_alpha, v_r, u, v, us, lhs, rhs, p

        if (mu < zero .or. ieee_is_nan(mu)) then
            k = -1
            return
        else if (mu == zero) then
            k = 0
            return
        end if

        if (mu < 30.0_dp) then
            p = exp(-mu)
            k = 0
            call random_number(u)
            do while (u > p)
                k = k + 1
                call random_number(v)
                u = u * v
            end do
            return
        end if

        ! PTRS transformed-rejection algorithm of Hoermann (1993).
        b = 0.931_dp + 2.53_dp * sqrt(mu)
        a = -0.059_dp + 0.02483_dp * b
        inv_alpha = 1.1239_dp + 1.1328_dp / (b - 3.4_dp)
        v_r = 0.9277_dp - 3.6224_dp / (b - 2.0_dp)
        do
            call random_number(u)
            u = u - 0.5_dp
            call random_number(v)
            us = 0.5_dp - abs(u)
            if (us <= zero) cycle
            k = floor((2.0_dp * a / us + b) * u + mu + 0.43_dp)
            if (us >= 0.07_dp .and. v <= v_r .and. k >= 0) return
            if (k < 0 .or. (us < 0.013_dp .and. v > us)) cycle
            lhs = log(v * inv_alpha / (a / (us * us) + b))
            rhs = -mu + real(k, dp) * log(mu) - log_gamma(real(k + 1, dp))
            if (lhs <= rhs) return
        end do
    end function poisson_random

end module delaporte_rng
