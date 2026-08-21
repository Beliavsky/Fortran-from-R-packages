! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_cvm
    use goftest_kinds, only : dp, pi
    use goftest_utils, only : sort_real, clamp01
    use goftest_special, only : bessel_k_frac
    implicit none
    private

    public :: cvm_cdf
    public :: cvm_quantile
    public :: cvm_statistic
    public :: cvm_test_uniform

contains

    real(dp) function ed2(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: z, b

        if (x <= 0.0_dp) then
            v = 0.0_dp
            return
        end if
        z = x * x / 4.0_dp
        b = bessel_k_frac(0.25_dp, z) + bessel_k_frac(0.75_dp, z)
        v = exp(-z) * b * sqrt(x**3 / (8.0_dp * pi))
    end function ed2

    real(dp) function ed3(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: z, k14, k34, k54, b

        if (x <= 0.0_dp) then
            v = 0.0_dp
            return
        end if
        z = x * x / 4.0_dp
        k14 = bessel_k_frac(0.25_dp, z)
        k34 = bessel_k_frac(0.75_dp, z)
        k54 = k34 + 0.5_dp * k14 / z
        b = 2.0_dp * k14 + 3.0_dp * k34 - k54
        v = exp(-z) * b * sqrt(x**5 / (32.0_dp * pi))
    end function ed3

    real(dp) function ak_over_factorial(k, x) result(v)
        integer, intent(in) :: k
        real(dp), intent(in) :: x
        real(dp) :: twosqrtx, fk1x, fk3x, fk5x, x34, x54, gf

        twosqrtx = 2.0_dp * sqrt(x)
        fk1x = real(4 * k + 1, dp) / twosqrtx
        fk3x = real(4 * k + 3, dp) / twosqrtx
        fk5x = real(4 * k + 5, dp) / twosqrtx
        x34 = x**0.75_dp
        x54 = x**1.25_dp
        if (k < 100) then
            gf = gamma(real(k, dp) + 0.5_dp) / gamma(real(k + 1, dp))
        else
            gf = exp(log_gamma(real(k, dp) + 0.5_dp) - log_gamma(real(k + 1, dp)))
        end if
        v = gf * (ed3(fk1x) / (72.0_dp * x54) + real(2 * k + 1, dp) * ( &
            ed2(fk3x) / (9.0_dp * x34) + real(2 * k + 3, dp) * ed3(fk5x) / (12.0_dp * x54) + &
            7.0_dp * (ed2(fk1x) + ed2(fk5x)) / (144.0_dp * x34)))
    end function ak_over_factorial

    real(dp) function cvm_asymptotic_cdf(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: total, q, term, coeff
        integer :: k

        if (x <= 0.0_dp) then
            v = 0.0_dp
            return
        end if
        total = 0.0_dp
        coeff = 1.0_dp
        do k = 0, 200
            if (k > 0) coeff = coeff * (real(k, dp) - 0.5_dp) / real(k, dp)
            q = real((4 * k + 1) * (4 * k + 1), dp) / (16.0_dp * x)
            term = coeff * sqrt(real(4 * k + 1, dp)) * exp(-q) * &
                bessel_k_frac(0.25_dp, q) / sqrt(x)
            total = total + term
            if (k > 10 .and. abs(term) < 1.0e-9_dp) exit
        end do
        v = clamp01(total / pi)
    end function cvm_asymptotic_cdf

    real(dp) function psi1(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: total, term
        integer :: k

        total = 0.0_dp
        do k = 0, 200
            term = -ak_over_factorial(k, x) / pi
            total = total + term
            if (k > 20 .and. abs(term) < 1.0e-9_dp) exit
        end do
        v = total + cvm_asymptotic_cdf(x) / 12.0_dp
    end function psi1

    real(dp) function cvm_cdf(q, n, lower_tail) result(p)
        real(dp), intent(in) :: q
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        real(dp) :: lower, upper
        logical :: lt
        integer :: nn

        lt = .true.
        if (present(lower_tail)) lt = lower_tail

        if (present(n)) then
            nn = min(100, max(1, n))
            lower = 1.0_dp / (12.0_dp * real(nn, dp))
            upper = real(nn, dp) / 3.0_dp
        else
            lower = 0.0_dp
            upper = huge(1.0_dp)
        end if

        if (q <= lower) then
            p = 0.0_dp
        else if (q >= upper) then
            p = 1.0_dp
        else if (present(n)) then
            p = cvm_asymptotic_cdf(q) + psi1(q) / real(n, dp)
        else
            p = cvm_asymptotic_cdf(q)
        end if
        p = clamp01(p)
        if (p < 2.0e-10_dp) p = 0.0_dp
        if (1.0_dp - p < 2.0e-10_dp) p = 1.0_dp
        if (.not. lt) p = 1.0_dp - p
    end function cvm_cdf

    real(dp) function cvm_quantile(prob, n, lower_tail) result(q)
        real(dp), intent(in) :: prob
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        real(dp) :: target, lo, hi, mid, fmid, lower, upper
        logical :: lt
        integer :: iter

        lt = .true.
        if (present(lower_tail)) lt = lower_tail
        target = prob
        if (.not. lt) target = 1.0_dp - target

        if (present(n)) then
            lower = 1.0_dp / (12.0_dp * real(n, dp))
            upper = real(n, dp) / 3.0_dp
        else
            lower = 0.0_dp
            upper = huge(1.0_dp)
        end if
        if (target <= 2.0e-10_dp) then
            q = lower
            return
        end if
        if (1.0_dp - target <= 2.0e-10_dp) then
            q = upper
            return
        end if
        lo = lower
        hi = 1.0_dp
        do
            if (present(n)) then
                fmid = cvm_cdf(hi, n=n)
            else
                fmid = cvm_cdf(hi)
            end if
            if (fmid >= target) exit
            hi = 2.0_dp * hi
            if (present(n)) hi = min(hi, upper)
            if (present(n) .and. hi >= upper) exit
        end do
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            if (present(n)) then
                fmid = cvm_cdf(mid, n=n)
            else
                fmid = cvm_cdf(mid)
            end if
            if (fmid < target) then
                lo = mid
            else
                hi = mid
            end if
        end do
        q = 0.5_dp * (lo + hi)
    end function cvm_quantile

    real(dp) function cvm_statistic(u) result(omega2)
        real(dp), intent(in) :: u(:)
        real(dp), allocatable :: x(:)
        integer :: i, n

        n = size(u)
        if (n == 0) then
            omega2 = 0.0_dp
            return
        end if
        allocate(x(n))
        x = u
        call sort_real(x)
        omega2 = 1.0_dp / (12.0_dp * real(n, dp))
        do i = 1, n
            omega2 = omega2 + (x(i) - real(2 * i - 1, dp) / (2.0_dp * real(n, dp)))**2
        end do
    end function cvm_statistic

    subroutine cvm_test_uniform(u, statistic, p_value)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: statistic, p_value
        integer :: n

        n = size(u)
        statistic = cvm_statistic(u)
        p_value = cvm_cdf(statistic, n=n, lower_tail=.false.)
    end subroutine cvm_test_uniform

end module goftest_cvm
