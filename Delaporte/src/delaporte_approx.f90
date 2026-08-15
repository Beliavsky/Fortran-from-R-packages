! Copyright (c) 2013, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! Monte-Carlo quantile approximation corresponding to qdelap(..., exact=FALSE)
! in the upstream R interface.

module delaporte_approx
    use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_positive_inf, &
        ieee_quiet_nan, ieee_value
    use delaporte_kinds, only : dp
    use delaporte_rng, only : gamma_random, poisson_random
    use delaporte_utils, only : zero, half, one
    implicit none
    private

    public :: qdelap_approx

contains

    subroutine qdelap_approx(p, alpha, beta, lambda, value, lower_tail, log_p, &
        nsim)
        real(dp), intent(in) :: p(:), alpha, beta, lambda
        real(dp), intent(out) :: value(size(p))
        logical, intent(in), optional :: lower_tail, log_p
        integer, intent(in), optional :: nsim
        integer, allocatable :: sample(:)
        integer :: i, nmc, expo
        real(dp) :: prob, mean_value
        logical :: lower, is_log

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        is_log = .false.
        if (present(log_p)) is_log = log_p

        if (alpha <= zero .or. beta <= zero .or. lambda <= zero) then
            value = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if

        if (present(nsim)) then
            nmc = max(1, nsim)
        else
            mean_value = alpha * beta + lambda
            expo = ceiling(log10(mean_value)) + 5
            if (expo >= 7) then
                nmc = 10000000
            else if (expo <= 0) then
                nmc = 1
            else
                nmc = 10 ** expo
            end if
        end if

        allocate(sample(nmc))
        do i = 1, nmc
            sample(i) = poisson_random(gamma_random(alpha, beta) + lambda)
        end do
        call sort_int(sample)

        do i = 1, size(p)
            prob = p(i)
            if (is_log) prob = exp(prob)
            if (.not. lower) prob = half - prob + half
            if (ieee_is_nan(prob) .or. prob < zero) then
                value(i) = ieee_value(prob, ieee_quiet_nan)
            else if (prob == zero) then
                value(i) = zero
            else if (prob >= one) then
                value(i) = ieee_value(prob, ieee_positive_inf)
            else
                value(i) = sample_quantile_type8(sample, prob)
            end if
        end do
    end subroutine qdelap_approx

    pure function sample_quantile_type8(x, p) result(q)
        integer, intent(in) :: x(:)
        real(dp), intent(in) :: p
        real(dp) :: q
        real(dp) :: h, gamma
        integer :: j, n

        n = size(x)
        h = (real(n, dp) + one / 3.0_dp) * p + one / 3.0_dp
        j = floor(h)
        gamma = h - real(j, dp)
        if (j <= 0) then
            q = real(x(1), dp)
        else if (j >= n) then
            q = real(x(n), dp)
        else
            q = (one - gamma) * real(x(j), dp) + gamma * real(x(j + 1), dp)
        end if
    end function sample_quantile_type8

    recursive subroutine quicksort_int(x, left, right)
        integer, intent(inout) :: x(:)
        integer, intent(in) :: left, right
        integer :: i, j, pivot, tmp

        i = left
        j = right
        pivot = x((left + right) / 2)
        do
            do while (x(i) < pivot)
                i = i + 1
            end do
            do while (x(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                tmp = x(i)
                x(i) = x(j)
                x(j) = tmp
                i = i + 1
                j = j - 1
            end if
            if (i > j) exit
        end do
        if (left < j) call quicksort_int(x, left, j)
        if (i < right) call quicksort_int(x, i, right)
    end subroutine quicksort_int

    subroutine sort_int(x)
        integer, intent(inout) :: x(:)

        if (size(x) > 1) call quicksort_int(x, 1, size(x))
    end subroutine sort_int

end module delaporte_approx
