! Copyright (c) 2016, Avraham Adler
! All rights reserved.
! SPDX-License-Identifier: BSD-2-Clause
!
! PMF/CDF/quantile algorithms adapted from src/delaporte.f90 in the
! Delaporte R package, with the R/C ABI removed and an idiomatic Fortran API.

module delaporte_distribution
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan, &
        ieee_positive_inf, ieee_quiet_nan, ieee_value
    use delaporte_kinds, only : dp, i64
    use delaporte_utils, only : zero, half, one, clamp01, log1p_delap, &
        recycle_index
    implicit none
    private

    real(dp), parameter :: max_exact_integer = real(huge(1_i64), dp)

    public :: ddelap, pdelap, qdelap
    public :: ddelap_vec, pdelap_vec, qdelap_vec

contains

    pure elemental function ddelap(x, alpha, beta, lambda, log_p) result(pmf)
        real(dp), intent(in) :: x, alpha, beta, lambda
        logical, intent(in), optional :: log_p
        real(dp) :: pmf
        real(dp) :: ii, kk, term
        integer(i64) :: i, k
        logical :: want_log

        want_log = .false.
        if (present(log_p)) want_log = log_p

        if (alpha <= zero .or. beta <= zero .or. lambda <= zero .or. &
            x < zero .or. ieee_is_nan(alpha + beta + lambda + x)) then
            pmf = ieee_value(x, ieee_quiet_nan)
            return
        end if

        if (.not. ieee_is_finite(x)) then
            pmf = zero
            if (want_log) pmf = log(pmf)
            return
        end if

        if (x >= max_exact_integer) then
            pmf = zero
            if (want_log) pmf = log(pmf)
            return
        end if

        k = floor(x, kind=i64)
        kk = real(k, dp)
        if (x /= kk) then
            pmf = zero
            if (want_log) pmf = log(pmf)
            return
        end if

        pmf = zero
        do i = 0_i64, k
            ii = real(i, dp)
            term = log_gamma(alpha + ii) + ii * log(beta) + &
                (kk - ii) * log(lambda) - lambda - log_gamma(alpha) - &
                log_gamma(ii + one) - (alpha + ii) * log1p_delap(beta) - &
                log_gamma(kk - ii + one)
            pmf = pmf + exp(term)
        end do
        pmf = clamp01(pmf)
        if (want_log) pmf = log(pmf)
    end function ddelap

    pure elemental function pdelap(q, alpha, beta, lambda, lower_tail, log_p) &
        result(cdf)
        real(dp), intent(in) :: q, alpha, beta, lambda
        logical, intent(in), optional :: lower_tail, log_p
        real(dp) :: cdf
        integer(i64) :: i, k
        logical :: lower, want_log

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        want_log = .false.
        if (present(log_p)) want_log = log_p

        if (alpha <= zero .or. beta <= zero .or. lambda <= zero .or. &
            q < zero .or. ieee_is_nan(alpha + beta + lambda + q)) then
            cdf = ieee_value(q, ieee_quiet_nan)
            return
        end if

        if (.not. ieee_is_finite(q) .or. q >= max_exact_integer) then
            cdf = one
        else
            k = floor(q, kind=i64)
            cdf = exp(-lambda) / ((beta + one) ** alpha)
            do i = 1_i64, k
                cdf = cdf + ddelap(real(i, dp), alpha, beta, lambda)
            end do
            cdf = clamp01(cdf)
        end if

        if (.not. lower) cdf = half - cdf + half
        if (want_log) cdf = log(cdf)
    end function pdelap

    pure elemental function qdelap(p, alpha, beta, lambda, lower_tail, log_p) &
        result(value)
        real(dp), intent(in) :: p, alpha, beta, lambda
        logical, intent(in), optional :: lower_tail, log_p
        real(dp) :: value
        real(dp) :: prob, testcdf
        logical :: lower, is_log

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        is_log = .false.
        if (present(log_p)) is_log = log_p

        prob = p
        if (is_log) prob = exp(prob)
        if (.not. lower) prob = half - prob + half

        if (alpha <= zero .or. beta <= zero .or. lambda <= zero .or. &
            prob < zero .or. ieee_is_nan(alpha + beta + lambda + prob)) then
            value = ieee_value(prob, ieee_quiet_nan)
        else if (prob >= one) then
            value = ieee_value(prob, ieee_positive_inf)
        else
            value = zero
            testcdf = exp(-lambda) / ((beta + one) ** alpha)
            do while (prob > testcdf)
                value = value + one
                testcdf = clamp01(testcdf + ddelap(value, alpha, beta, lambda))
            end do
        end if
    end function qdelap

    subroutine ddelap_vec(x, alpha, beta, lambda, pmf, log_p)
        real(dp), intent(in) :: x(:), alpha(:), beta(:), lambda(:)
        real(dp), intent(out) :: pmf(size(x))
        logical, intent(in), optional :: log_p
        integer :: i
        logical :: want_log

        want_log = .false.
        if (present(log_p)) want_log = log_p
        do i = 1, size(x)
            pmf(i) = ddelap(x(i), alpha(recycle_index(i, size(alpha))), &
                beta(recycle_index(i, size(beta))), &
                lambda(recycle_index(i, size(lambda))), want_log)
        end do
    end subroutine ddelap_vec

    subroutine pdelap_vec(q, alpha, beta, lambda, cdf, lower_tail, log_p)
        real(dp), intent(in) :: q(:), alpha(:), beta(:), lambda(:)
        real(dp), intent(out) :: cdf(size(q))
        logical, intent(in), optional :: lower_tail, log_p
        integer :: i
        logical :: lower, want_log

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        want_log = .false.
        if (present(log_p)) want_log = log_p
        do i = 1, size(q)
            cdf(i) = pdelap(q(i), alpha(recycle_index(i, size(alpha))), &
                beta(recycle_index(i, size(beta))), &
                lambda(recycle_index(i, size(lambda))), lower, want_log)
        end do
    end subroutine pdelap_vec

    subroutine qdelap_vec(p, alpha, beta, lambda, value, lower_tail, log_p)
        real(dp), intent(in) :: p(:), alpha(:), beta(:), lambda(:)
        real(dp), intent(out) :: value(size(p))
        logical, intent(in), optional :: lower_tail, log_p
        integer :: i
        logical :: lower, is_log

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        is_log = .false.
        if (present(log_p)) is_log = log_p
        do i = 1, size(p)
            value(i) = qdelap(p(i), alpha(recycle_index(i, size(alpha))), &
                beta(recycle_index(i, size(beta))), &
                lambda(recycle_index(i, size(lambda))), lower, is_log)
        end do
    end subroutine qdelap_vec

end module delaporte_distribution
