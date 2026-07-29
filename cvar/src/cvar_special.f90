! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of cvar 0.6 by Georgi N. Boshnakov.
module cvar_special
    use cvar_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
    real(dp), parameter :: log_sqrt_two_pi = 0.5_dp * log(2.0_dp * pi)

    public :: normal_pdf, normal_cdf, normal_quantile
    public :: student_t_pdf, student_t_cdf, student_t_quantile
    public :: std_student_t_pdf, std_student_t_cdf, std_student_t_quantile
    public :: ged_pdf, ged_cdf, ged_quantile, ged_scale
    public :: regularized_gamma_p, regularized_beta

contains

    pure function normal_pdf(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value
        value = exp(-0.5_dp * x * x - log_sqrt_two_pi)
    end function normal_pdf

    pure function normal_cdf(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value
        value = 0.5_dp * erfc(-x / sqrt_two)
    end function normal_cdf

    pure function normal_quantile(p) result(value)
        real(dp), intent(in) :: p
        real(dp) :: value
        real(dp) :: q, r
        real(dp), parameter :: a(6) = [ &
            -3.969683028665376d1, 2.209460984245205d2, &
            -2.759285104469687d2, 1.383577518672690d2, &
            -3.066479806614716d1, 2.506628277459239d0 ]
        real(dp), parameter :: b(5) = [ &
            -5.447609879822406d1, 1.615858368580409d2, &
            -1.556989798598866d2, 6.680131188771972d1, &
            -1.328068155288572d1 ]
        real(dp), parameter :: c(6) = [ &
            -7.784894002430293d-3, -3.223964580411365d-1, &
            -2.400758277161838d0, -2.549732539343734d0, &
             4.374664141464968d0, 2.938163982698783d0 ]
        real(dp), parameter :: d(4) = [ &
             7.784695709041462d-3, 3.224671290700398d-1, &
             2.445134137142996d0, 3.754408661907416d0 ]
        real(dp), parameter :: p_low = 0.02425_dp
        real(dp), parameter :: p_high = 1.0_dp - p_low

        if (p <= 0.0_dp) then
            value = -huge(1.0_dp)
        else if (p >= 1.0_dp) then
            value = huge(1.0_dp)
        else if (p < p_low) then
            q = sqrt(-2.0_dp * log(p))
            value = (((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
                    ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
        else if (p <= p_high) then
            q = p - 0.5_dp
            r = q * q
            value = (((((a(1) * r + a(2)) * r + a(3)) * r + a(4)) * r + a(5)) * r + a(6)) * q / &
                    (((((b(1) * r + b(2)) * r + b(3)) * r + b(4)) * r + b(5)) * r + 1.0_dp)
        else
            q = sqrt(-2.0_dp * log(1.0_dp - p))
            value = -(((((c(1) * q + c(2)) * q + c(3)) * q + c(4)) * q + c(5)) * q + c(6)) / &
                     ((((d(1) * q + d(2)) * q + d(3)) * q + d(4)) * q + 1.0_dp)
        end if

        if (p > 0.0_dp .and. p < 1.0_dp) then
            value = value - (normal_cdf(value) - p) / max(normal_pdf(value), tiny(1.0_dp))
        end if
    end function normal_quantile

    pure function regularized_beta(a, b, x) result(value)
        real(dp), intent(in) :: a, b, x
        real(dp) :: value, bt

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        else if (x >= 1.0_dp) then
            value = 1.0_dp
            return
        end if

        bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + &
                 a * log(x) + b * log(1.0_dp - x))
        if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
            value = bt * beta_cf(a, b, x) / a
        else
            value = 1.0_dp - bt * beta_cf(b, a, 1.0_dp - x) / b
        end if
        value = min(1.0_dp, max(0.0_dp, value))
    end function regularized_beta

    pure function beta_cf(a, b, x) result(value)
        real(dp), intent(in) :: a, b, x
        real(dp) :: value
        integer :: m, m2
        real(dp) :: aa, c, d, del, h, qab, qam, qap
        real(dp), parameter :: eps = 3.0e-14_dp
        real(dp), parameter :: fpmin = 1.0e-300_dp

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab * x / qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp / d
        h = d

        do m = 1, 400
            m2 = 2 * m
            aa = real(m, dp) * (b - real(m, dp)) * x / &
                 ((qam + real(m2, dp)) * (a + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            h = h * d * c

            aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x / &
                 ((a + real(m2, dp)) * (qap + real(m2, dp)))
            d = 1.0_dp + aa * d
            if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        value = h
    end function beta_cf

    pure function regularized_gamma_p(a, x) result(value)
        real(dp), intent(in) :: a, x
        real(dp) :: value
        integer :: n
        real(dp) :: ap, del, sum, b, c, d, h, an
        real(dp), parameter :: eps = 3.0e-14_dp
        real(dp), parameter :: fpmin = 1.0e-300_dp

        if (x <= 0.0_dp) then
            value = 0.0_dp
            return
        end if
        if (a <= 0.0_dp) then
            value = 0.0_dp
            return
        end if

        if (x < a + 1.0_dp) then
            ap = a
            sum = 1.0_dp / a
            del = sum
            do n = 1, 1000
                ap = ap + 1.0_dp
                del = del * x / ap
                sum = sum + del
                if (abs(del) <= abs(sum) * eps) exit
            end do
            value = sum * exp(-x + a * log(x) - log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp / fpmin
            d = 1.0_dp / b
            h = d
            do n = 1, 1000
                an = -real(n, dp) * (real(n, dp) - a)
                b = b + 2.0_dp
                d = an * d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an / c
                if (abs(c) < fpmin) c = fpmin
                d = 1.0_dp / d
                del = d * c
                h = h * del
                if (abs(del - 1.0_dp) <= eps) exit
            end do
            value = 1.0_dp - exp(-x + a * log(x) - log_gamma(a)) * h
        end if
        value = min(1.0_dp, max(0.0_dp, value))
    end function regularized_gamma_p

    pure function student_t_pdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value
        value = exp(log_gamma(0.5_dp * (nu + 1.0_dp)) - log_gamma(0.5_dp * nu) - &
                    0.5_dp * log(nu * pi) - 0.5_dp * (nu + 1.0_dp) * log(1.0_dp + x * x / nu))
    end function student_t_pdf

    pure function student_t_cdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value, ib, z
        if (abs(x) <= tiny(1.0_dp)) then
            value = 0.5_dp
            return
        end if
        z = nu / (nu + x * x)
        ib = regularized_beta(0.5_dp * nu, 0.5_dp, z)
        if (x > 0.0_dp) then
            value = 1.0_dp - 0.5_dp * ib
        else
            value = 0.5_dp * ib
        end if
    end function student_t_cdf

    pure function student_t_quantile(p, nu) result(value)
        real(dp), intent(in) :: p, nu
        real(dp) :: value, lo, hi, mid
        integer :: iter

        if (p <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        else if (p >= 1.0_dp) then
            value = huge(1.0_dp)
            return
        else if (abs(p - 0.5_dp) <= epsilon(1.0_dp)) then
            value = 0.0_dp
            return
        end if

        lo = -1.0_dp
        hi = 1.0_dp
        do while (student_t_cdf(lo, nu) > p)
            lo = 2.0_dp * lo
            if (lo < -1.0e12_dp) exit
        end do
        do while (student_t_cdf(hi, nu) < p)
            hi = 2.0_dp * hi
            if (hi > 1.0e12_dp) exit
        end do
        do iter = 1, 200
            mid = 0.5_dp * (lo + hi)
            if (student_t_cdf(mid, nu) < p) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 2.0e-13_dp * (1.0_dp + abs(mid))) exit
        end do
        value = 0.5_dp * (lo + hi)
    end function student_t_quantile

    pure function std_student_t_pdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value, scale
        scale = sqrt((nu - 2.0_dp) / nu)
        value = student_t_pdf(x / scale, nu) / scale
    end function std_student_t_pdf

    pure function std_student_t_cdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value, scale
        scale = sqrt((nu - 2.0_dp) / nu)
        value = student_t_cdf(x / scale, nu)
    end function std_student_t_cdf

    pure function std_student_t_quantile(p, nu) result(value)
        real(dp), intent(in) :: p, nu
        real(dp) :: value
        value = sqrt((nu - 2.0_dp) / nu) * student_t_quantile(p, nu)
    end function std_student_t_quantile

    pure function ged_scale(nu) result(value)
        real(dp), intent(in) :: nu
        real(dp) :: value
        value = sqrt(2.0_dp**(-2.0_dp / nu) * exp(log_gamma(1.0_dp / nu) - &
                     log_gamma(3.0_dp / nu)))
    end function ged_scale

    pure function ged_pdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value, lambda
        lambda = ged_scale(nu)
        value = exp(log(nu) - log(lambda) - (1.0_dp + 1.0_dp / nu) * log(2.0_dp) - &
                    log_gamma(1.0_dp / nu) - 0.5_dp * abs(x / lambda)**nu)
    end function ged_pdf

    pure function ged_cdf(x, nu) result(value)
        real(dp), intent(in) :: x, nu
        real(dp) :: value, p, z, lambda
        if (abs(x) <= tiny(1.0_dp)) then
            value = 0.5_dp
            return
        end if
        lambda = ged_scale(nu)
        z = 0.5_dp * abs(x / lambda)**nu
        p = regularized_gamma_p(1.0_dp / nu, z)
        if (x > 0.0_dp) then
            value = 0.5_dp * (1.0_dp + p)
        else
            value = 0.5_dp * (1.0_dp - p)
        end if
    end function ged_cdf

    pure function ged_quantile(p, nu) result(value)
        real(dp), intent(in) :: p, nu
        real(dp) :: value, lo, hi, mid
        integer :: iter
        if (p <= 0.0_dp) then
            value = -huge(1.0_dp)
            return
        else if (p >= 1.0_dp) then
            value = huge(1.0_dp)
            return
        else if (abs(p - 0.5_dp) <= epsilon(1.0_dp)) then
            value = 0.0_dp
            return
        end if
        lo = -1.0_dp
        hi = 1.0_dp
        do while (ged_cdf(lo, nu) > p)
            lo = 2.0_dp * lo
            if (lo < -1.0e12_dp) exit
        end do
        do while (ged_cdf(hi, nu) < p)
            hi = 2.0_dp * hi
            if (hi > 1.0e12_dp) exit
        end do
        do iter = 1, 200
            mid = 0.5_dp * (lo + hi)
            if (ged_cdf(mid, nu) < p) then
                lo = mid
            else
                hi = mid
            end if
            if (abs(hi - lo) <= 2.0e-13_dp * (1.0_dp + abs(mid))) exit
        end do
        value = 0.5_dp * (lo + hi)
    end function ged_quantile

end module cvar_special
