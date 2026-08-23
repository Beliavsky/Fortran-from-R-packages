module skewt_special
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_negative_inf, ieee_is_finite
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: max_iter = 400
    real(dp), parameter :: eps = 4.0_dp * epsilon(1.0_dp)
    real(dp), parameter :: fpmin = tiny(1.0_dp) / eps

    public :: student_t_pdf, student_t_cdf, student_t_quantile
    public :: student_t_random

contains

    pure real(dp) function nan_value() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_value

    pure real(dp) function student_t_pdf(x, df) result(f)
        real(dp), intent(in) :: x, df
        real(dp) :: logf

        if (df <= 0.0_dp) then
            f = nan_value()
            return
        end if

        logf = log_gamma(0.5_dp * (df + 1.0_dp)) - log_gamma(0.5_dp * df) &
            - 0.5_dp * (log(df) + log(pi)) &
            - 0.5_dp * (df + 1.0_dp) * log(1.0_dp + x * x / df)
        f = exp(logf)
    end function student_t_pdf

    real(dp) function student_t_cdf(x, df) result(p)
        real(dp), intent(in) :: x, df
        real(dp) :: z, ib

        if (df <= 0.0_dp) then
            p = nan_value()
            return
        end if
        if (.not. ieee_is_finite(x)) then
            if (x < 0.0_dp) then
                p = 0.0_dp
            else
                p = 1.0_dp
            end if
            return
        end if

        z = df / (df + x * x)
        ib = regularized_beta(z, 0.5_dp * df, 0.5_dp)
        if (x > 0.0_dp) then
            p = 1.0_dp - 0.5_dp * ib
        else
            p = 0.5_dp * ib
        end if
        p = min(1.0_dp, max(0.0_dp, p))
    end function student_t_cdf

    recursive real(dp) function student_t_quantile(p, df) result(x)
        real(dp), intent(in) :: p, df
        real(dp) :: lo, hi, mid, pmid
        integer :: iter

        if (df <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = nan_value()
            return
        end if
        if (p <= 0.0_dp) then
            x = ieee_value(x, ieee_negative_inf)
            return
        else if (p >= 1.0_dp) then
            x = ieee_value(x, ieee_positive_inf)
            return
        else if (p <= 0.5_dp .and. p >= 0.5_dp) then
            x = 0.0_dp
            return
        end if

        if (p < 0.5_dp) then
            x = -student_t_quantile(1.0_dp - p, df)
            return
        end if

        lo = 0.0_dp
        hi = 1.0_dp
        do while (student_t_cdf(hi, df) < p)
            hi = 2.0_dp * hi
            if (hi > sqrt(huge(1.0_dp))) then
                x = ieee_value(x, ieee_positive_inf)
                return
            end if
        end do

        do iter = 1, 180
            mid = 0.5_dp * (lo + hi)
            pmid = student_t_cdf(mid, df)
            if (pmid < p) then
                lo = mid
            else
                hi = mid
            end if
            if (hi - lo <= 8.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(mid))) exit
        end do
        x = 0.5_dp * (lo + hi)
    end function student_t_quantile

    real(dp) function student_t_random(df) result(x)
        real(dp), intent(in) :: df
        real(dp) :: z, g

        if (df <= 0.0_dp) then
            x = nan_value()
            return
        end if
        z = normal_random()
        g = gamma_random(0.5_dp * df)
        x = z / sqrt((2.0_dp * g) / df)
    end function student_t_random

    real(dp) function regularized_beta(x, a, b) result(val)
        real(dp), intent(in) :: x, a, b
        real(dp) :: bt

        if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
            val = nan_value()
            return
        end if
        if (x <= 0.0_dp) then
            val = 0.0_dp
            return
        else if (x >= 1.0_dp) then
            val = 1.0_dp
            return
        end if

        bt = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) &
            + a * log(x) + b * log(1.0_dp - x))
        if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
            val = bt * beta_cf(a, b, x) / a
        else
            val = 1.0_dp - bt * beta_cf(b, a, 1.0_dp - x) / b
        end if
        val = min(1.0_dp, max(0.0_dp, val))
    end function regularized_beta

    real(dp) function beta_cf(a, b, x) result(h)
        real(dp), intent(in) :: a, b, x
        real(dp) :: qab, qap, qam, c, d, aa, del
        integer :: m, m2

        qab = a + b
        qap = a + 1.0_dp
        qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab * x / qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp / d
        h = d

        do m = 1, max_iter
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
    end function beta_cf

    real(dp) function normal_random() result(z)
        real(dp) :: u1, u2
        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * pi * u2)
    end function normal_random

    recursive real(dp) function gamma_random(shape) result(g)
        real(dp), intent(in) :: shape
        real(dp) :: d, c, x, v, u

        if (shape <= 0.0_dp) then
            g = nan_value()
            return
        end if
        if (shape < 1.0_dp) then
            call random_number(u)
            u = max(u, tiny(1.0_dp))
            g = gamma_random(shape + 1.0_dp) * u**(1.0_dp / shape)
            return
        end if

        d = shape - 1.0_dp / 3.0_dp
        c = 1.0_dp / sqrt(9.0_dp * d)
        do
            do
                x = normal_random()
                v = 1.0_dp + c * x
                if (v > 0.0_dp) exit
            end do
            v = v * v * v
            call random_number(u)
            if (u < 1.0_dp - 0.0331_dp * x**4) exit
            if (log(max(u, tiny(1.0_dp))) < 0.5_dp * x * x &
                + d * (1.0_dp - v + log(v))) exit
        end do
        g = d * v
    end function gamma_random

end module skewt_special
