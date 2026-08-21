module circstats_special
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use circstats_kinds, only: dp, pi
    implicit none
    private
    public :: normal_pdf, normal_cdf, gamma_p, gamma_q, chi_square_cdf
    public :: chi_square_quantile, i0, i1, ip, i0e, i1e, log_i0

contains

    pure elemental real(dp) function normal_pdf(x) result(y)
        real(dp), intent(in) :: x
        y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
    end function normal_pdf

    pure elemental real(dp) function normal_cdf(x) result(y)
        real(dp), intent(in) :: x
        y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
    end function normal_cdf

    pure elemental real(dp) function i0e(x) result(ans)
        real(dp), intent(in) :: x
        real(dp) :: ax, y, z
        ax = abs(x)
        if (ax < 3.75_dp) then
            y = (ax/3.75_dp)**2
            ans = 1.0_dp + y*(3.5156229_dp + y*(3.0899424_dp + y*(1.2067492_dp + &
                y*(0.2659732_dp + y*(0.0360768_dp + y*0.0045813_dp)))))
            ans = ans*exp(-ax)
        else
            z = 3.75_dp/ax
            ans = 0.39894228_dp + z*(0.01328592_dp + z*(0.00225319_dp + z*(-0.00157565_dp + &
                z*(0.00916281_dp + z*(-0.02057706_dp + z*(0.02635537_dp + &
                z*(-0.01647633_dp + z*0.00392377_dp)))))))
            ans = ans/sqrt(ax)
        end if
    end function i0e

    pure elemental real(dp) function i1e(x) result(ans)
        real(dp), intent(in) :: x
        real(dp) :: ax, y, z
        ax = abs(x)
        if (ax < 3.75_dp) then
            y = (ax/3.75_dp)**2
            ans = ax*(0.5_dp + y*(0.87890594_dp + y*(0.51498869_dp + y*(0.15084934_dp + &
                y*(0.02658733_dp + y*(0.00301532_dp + y*0.00032411_dp))))))
            ans = ans*exp(-ax)
        else
            z = 3.75_dp/ax
            ans = 0.39894228_dp + z*(-0.03988024_dp + z*(-0.00362018_dp + z*(0.00163801_dp + &
                z*(-0.01031555_dp + z*(0.02282967_dp + z*(-0.02895312_dp + &
                z*(0.01787654_dp - z*0.00420059_dp)))))))
            ans = ans/sqrt(ax)
        end if
        if (x < 0.0_dp) ans = -ans
    end function i1e

    pure elemental real(dp) function i0(x) result(ans)
        real(dp), intent(in) :: x
        if (abs(x) > log(huge(1.0_dp))) then
            ans = ieee_value(ans, ieee_positive_inf)
        else
            ans = i0e(x)*exp(abs(x))
        end if
    end function i0

    pure elemental real(dp) function i1(x) result(ans)
        real(dp), intent(in) :: x
        if (abs(x) > log(huge(1.0_dp))) then
            ans = sign(ieee_value(ans, ieee_positive_inf), x)
        else
            ans = i1e(x)*exp(abs(x))
        end if
    end function i1

    pure real(dp) function ip(n, x) result(ans)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: scaled, term, sumv, xx
        integer :: m
        if (n < 0) then
            ans = 0.0_dp
            return
        end if
        if (n == 0) then
            ans = i0(x)
            return
        end if
        if (n == 1) then
            ans = i1(x)
            return
        end if
        xx = abs(x)
        if (xx <= tiny(1.0_dp)) then
            ans = 0.0_dp
            return
        end if
        term = exp(real(n,dp)*log(0.5_dp*xx) - log_gamma(real(n+1,dp)) - xx)
        sumv = term
        do m = 1, 200000
            term = term*(0.25_dp*xx*xx)/(real(m,dp)*real(n+m,dp))
            sumv = sumv + term
            if (abs(term) <= 2.0_dp*epsilon(1.0_dp)*max(abs(sumv),tiny(1.0_dp))) exit
        end do
        scaled = sumv
        if (xx > log(huge(1.0_dp))) then
            ans = ieee_value(ans, ieee_positive_inf)
        else
            ans = scaled*exp(xx)
        end if
        if (x < 0.0_dp .and. mod(n,2) == 1) ans = -ans
    end function ip

    pure elemental real(dp) function log_i0(x) result(ans)
        real(dp), intent(in) :: x
        ans = log(i0e(x)) + abs(x)
    end function log_i0

    pure real(dp) function gamma_p(a, x) result(p)
        real(dp), intent(in) :: a, x
        real(dp), parameter :: eps = 3.0e-14_dp
        real(dp), parameter :: fpmin = tiny(1.0_dp)/eps
        integer, parameter :: itmax = 10000
        real(dp) :: ap, del, sumv, b, c, d, h, an
        integer :: n
        if (a <= 0.0_dp .or. x < 0.0_dp) then
            p = 0.0_dp
            return
        end if
        if (x <= tiny(1.0_dp)) then
            p = 0.0_dp
            return
        end if
        if (x < a + 1.0_dp) then
            ap = a
            sumv = 1.0_dp/a
            del = sumv
            do n = 1, itmax
                ap = ap + 1.0_dp
                del = del*x/ap
                sumv = sumv + del
                if (abs(del) < abs(sumv)*eps) exit
            end do
            p = sumv*exp(-x + a*log(x) - log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp/fpmin
            d = 1.0_dp/b
            h = d
            do n = 1, itmax
                an = -real(n,dp)*(real(n,dp)-a)
                b = b + 2.0_dp
                d = an*d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an/c
                if (abs(c) < fpmin) c = fpmin
                d = 1.0_dp/d
                del = d*c
                h = h*del
                if (abs(del - 1.0_dp) < eps) exit
            end do
            p = 1.0_dp - exp(-x + a*log(x) - log_gamma(a))*h
        end if
        p = max(0.0_dp, min(1.0_dp, p))
    end function gamma_p

    pure real(dp) function gamma_q(a, x) result(q)
        real(dp), intent(in) :: a, x
        q = 1.0_dp - gamma_p(a, x)
    end function gamma_q

    pure real(dp) function chi_square_cdf(x, df) result(p)
        real(dp), intent(in) :: x
        integer, intent(in) :: df
        if (x <= 0.0_dp) then
            p = 0.0_dp
        else
            p = gamma_p(0.5_dp*real(df,dp), 0.5_dp*x)
        end if
    end function chi_square_cdf

    pure real(dp) function chi_square_quantile(prob, df) result(x)
        real(dp), intent(in) :: prob
        integer, intent(in) :: df
        real(dp) :: lo, hi, mid
        integer :: iter
        if (prob <= 0.0_dp) then
            x = 0.0_dp
            return
        end if
        if (prob >= 1.0_dp) then
            x = huge(1.0_dp)
            return
        end if
        lo = 0.0_dp
        hi = max(1.0_dp, real(df,dp))
        do while (chi_square_cdf(hi,df) < prob)
            hi = 2.0_dp*hi
            if (hi > 0.25_dp*huge(1.0_dp)) exit
        end do
        do iter = 1, 200
            mid = 0.5_dp*(lo+hi)
            if (chi_square_cdf(mid,df) < prob) then
                lo = mid
            else
                hi = mid
            end if
            if (hi-lo <= 1.0e-13_dp*max(1.0_dp,mid)) exit
        end do
        x = 0.5_dp*(lo+hi)
    end function chi_square_quantile
end module circstats_special
