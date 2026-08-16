! SPDX-License-Identifier: LGPL-3.0-or-later
! Standalone numerical support for the PDQutils translation.
module pdqutils_special
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use pdqutils_kinds, only : dp, pi, sqrt2
    implicit none
    private

    public :: nan_dp, pos_inf_dp, neg_inf_dp
    public :: factorial_dp, binomial_dp
    public :: normal_pdf, normal_cdf, normal_quantile
    public :: gamma_pdf, gamma_cdf, beta_pdf, beta_cdf
    public :: log_beta_dp

contains

    pure real(dp) function nan_dp() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_dp

    pure real(dp) function pos_inf_dp() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function pos_inf_dp

    pure real(dp) function neg_inf_dp() result(x)
        x = ieee_value(0.0_dp, ieee_negative_inf)
    end function neg_inf_dp

    pure real(dp) function factorial_dp(n) result(v)
        integer, intent(in) :: n
        integer :: j
        if (n < 0) then
            v = nan_dp()
            return
        end if
        v = 1.0_dp
        do j = 2, n
            v = v*real(j,dp)
        end do
    end function factorial_dp

    pure real(dp) function binomial_dp(n,k) result(v)
        integer, intent(in) :: n, k
        integer :: j, kk
        if (k < 0 .or. k > n) then
            v = 0.0_dp
            return
        end if
        kk = min(k,n-k)
        v = 1.0_dp
        do j = 1, kk
            v = v*real(n-kk+j,dp)/real(j,dp)
        end do
    end function binomial_dp

    pure real(dp) function normal_pdf(x) result(v)
        real(dp), intent(in) :: x
        v = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
    end function normal_pdf

    pure real(dp) function normal_cdf(x, lower_tail) result(v)
        real(dp), intent(in) :: x
        logical, intent(in), optional :: lower_tail
        logical :: lower
        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        if (lower) then
            v = 0.5_dp*erfc(-x/sqrt2)
        else
            v = 0.5_dp*erfc(x/sqrt2)
        end if
    end function normal_cdf

    pure real(dp) function expm1_local(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: term
        integer :: j
        if (abs(x) > 1.0e-5_dp) then
            v = exp(x)-1.0_dp
        else
            v = x
            term = x
            do j = 2, 24
                term = term*x/real(j,dp)
                v = v + term
            end do
        end if
    end function expm1_local

    pure real(dp) function log1p_local(x) result(v)
        real(dp), intent(in) :: x
        real(dp) :: term, add
        integer :: j
        if (x <= -1.0_dp) then
            if (x == -1.0_dp) then
                v = neg_inf_dp()
            else
                v = nan_dp()
            end if
        else if (abs(x) > 1.0e-5_dp) then
            v = log(1.0_dp+x)
        else
            v = 0.0_dp
            term = x
            do j = 1, 48
                add = term/real(j,dp)
                if (mod(j,2) == 0) add = -add
                v = v + add
                term = term*x
            end do
        end if
    end function log1p_local

    pure real(dp) function normal_quantile(p, lower_tail, log_p) result(x)
        ! Acklam inverse-normal approximation with a Halley correction.
        real(dp), intent(in) :: p
        logical, intent(in), optional :: lower_tail, log_p
        logical :: lower, lp
        real(dp) :: pp, q, r, e, u, sign_out, log_small
        real(dp), parameter :: a1=-3.969683028665376e1_dp, a2=2.209460984245205e2_dp
        real(dp), parameter :: a3=-2.759285104469687e2_dp, a4=1.383577518672690e2_dp
        real(dp), parameter :: a5=-3.066479806614716e1_dp, a6=2.506628277459239_dp
        real(dp), parameter :: b1=-5.447609879822406e1_dp, b2=1.615858368580409e2_dp
        real(dp), parameter :: b3=-1.556989798598866e2_dp, b4=6.680131188771972e1_dp
        real(dp), parameter :: b5=-1.328068155288572e1_dp
        real(dp), parameter :: c1=-7.784894002430293e-3_dp, c2=-3.223964580411365e-1_dp
        real(dp), parameter :: c3=-2.400758277161838_dp, c4=-2.549732539343734_dp
        real(dp), parameter :: c5=4.374664141464968_dp, c6=2.938163982698783_dp
        real(dp), parameter :: d1=7.784695709041462e-3_dp, d2=3.224671290700398e-1_dp
        real(dp), parameter :: d3=2.445134137142996_dp, d4=3.754408661907416_dp
        real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow

        lower = .true.; lp = .false.
        if (present(lower_tail)) lower = lower_tail
        if (present(log_p)) lp = log_p
        sign_out = merge(1.0_dp,-1.0_dp,lower)

        if (lp) then
            if (p > 0.0_dp) then
                x = nan_dp(); return
            else if (p == 0.0_dp) then
                x = sign_out*pos_inf_dp(); return
            else if (p == neg_inf_dp()) then
                x = sign_out*neg_inf_dp(); return
            end if
            if (p <= log(0.5_dp)) then
                if (p < log(plow)) then
                    q = sqrt(-2.0_dp*p)
                    x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                        ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
                    x = sign_out*x
                    return
                end if
                pp = exp(p)
            else
                log_small = log(-expm1_local(p))
                if (log_small < log(plow)) then
                    q = sqrt(-2.0_dp*log_small)
                    x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                        ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
                    x = sign_out*x
                    return
                end if
                pp = exp(p)
            end if
        else
            pp = p
            if (pp < 0.0_dp .or. pp > 1.0_dp) then
                x = nan_dp(); return
            else if (pp == 0.0_dp) then
                x = sign_out*neg_inf_dp(); return
            else if (pp == 1.0_dp) then
                x = sign_out*pos_inf_dp(); return
            end if
        end if

        if (pp < plow) then
            q = sqrt(-2.0_dp*log(pp))
            x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        else if (pp <= phigh) then
            q = pp-0.5_dp
            r = q*q
            x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
                (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
        else
            q = sqrt(-2.0_dp*log1p_local(-pp))
            x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        end if
        if (abs(x) < 8.0_dp) then
            e = normal_cdf(x)-pp
            u = e/normal_pdf(x)
            x = x-u/(1.0_dp+0.5_dp*x*u)
        end if
        x = sign_out*x
    end function normal_quantile

    pure real(dp) function gamma_pdf(x, shape) result(v)
        real(dp), intent(in) :: x, shape
        if (shape <= 0.0_dp .or. x < 0.0_dp) then
            v = nan_dp()
        else if (x == 0.0_dp) then
            if (shape < 1.0_dp) then
                v = pos_inf_dp()
            else if (shape == 1.0_dp) then
                v = 1.0_dp
            else
                v = 0.0_dp
            end if
        else
            v = exp((shape-1.0_dp)*log(x)-x-log_gamma(shape))
        end if
    end function gamma_pdf

    pure real(dp) function regularized_gamma_p(a,x) result(p)
        real(dp), intent(in) :: a, x
        integer, parameter :: itmax=10000
        real(dp), parameter :: eps=4.0_dp*epsilon(1.0_dp), fpmin=tiny(1.0_dp)/eps
        real(dp) :: ap, del, sumv, b, c, d, h, an, gln, q
        integer :: n
        if (a <= 0.0_dp .or. x < 0.0_dp) then
            p = nan_dp(); return
        end if
        if (x == 0.0_dp) then
            p = 0.0_dp; return
        end if
        gln = log_gamma(a)
        if (x < a+1.0_dp) then
            ap = a
            sumv = 1.0_dp/a
            del = sumv
            do n=1,itmax
                ap = ap+1.0_dp
                del = del*x/ap
                sumv = sumv+del
                if (abs(del) <= abs(sumv)*eps) exit
            end do
            p = sumv*exp(-x+a*log(x)-gln)
        else
            b = x+1.0_dp-a
            c = 1.0_dp/fpmin
            d = 1.0_dp/b
            h = d
            do n=1,itmax
                an = -real(n,dp)*(real(n,dp)-a)
                b = b+2.0_dp
                d = an*d+b
                if (abs(d) < fpmin) d=fpmin
                c = b+an/c
                if (abs(c) < fpmin) c=fpmin
                d = 1.0_dp/d
                del = d*c
                h = h*del
                if (abs(del-1.0_dp) <= eps) exit
            end do
            q = exp(-x+a*log(x)-gln)*h
            p = 1.0_dp-q
        end if
        p = min(1.0_dp,max(0.0_dp,p))
    end function regularized_gamma_p

    pure real(dp) function gamma_cdf(x, shape) result(v)
        real(dp), intent(in) :: x, shape
        if (shape <= 0.0_dp) then
            v = nan_dp()
        else if (x <= 0.0_dp) then
            v = 0.0_dp
        else
            v = regularized_gamma_p(shape,x)
        end if
    end function gamma_cdf

    pure real(dp) function log_beta_dp(a,b) result(v)
        real(dp), intent(in) :: a,b
        if (a <= 0.0_dp .or. b <= 0.0_dp) then
            v = nan_dp()
        else
            v = log_gamma(a)+log_gamma(b)-log_gamma(a+b)
        end if
    end function log_beta_dp

    pure real(dp) function beta_pdf(x,a,b) result(v)
        real(dp), intent(in) :: x,a,b
        if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
            v = nan_dp()
        else if (x == 0.0_dp) then
            if (a < 1.0_dp) then
                v = pos_inf_dp()
            else if (a == 1.0_dp) then
                v = b
            else
                v = 0.0_dp
            end if
        else if (x == 1.0_dp) then
            if (b < 1.0_dp) then
                v = pos_inf_dp()
            else if (b == 1.0_dp) then
                v = a
            else
                v = 0.0_dp
            end if
        else
            v = exp((a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x)-log_beta_dp(a,b))
        end if
    end function beta_pdf

    pure real(dp) function betacf(a,b,x) result(h)
        real(dp), intent(in) :: a,b,x
        integer, parameter :: maxit=10000
        real(dp), parameter :: eps=4.0_dp*epsilon(1.0_dp), fpmin=tiny(1.0_dp)/eps
        real(dp) :: qab,qap,qam,c,d,del,aa
        integer :: m,m2
        qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
        c=1.0_dp
        d=1.0_dp-qab*x/qap
        if(abs(d)<fpmin)d=fpmin
        d=1.0_dp/d
        h=d
        do m=1,maxit
            m2=2*m
            aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
            d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
            c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d; h=h*d*c
            aa=-(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
            d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
            c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
            d=1.0_dp/d; del=d*c; h=h*del
            if(abs(del-1.0_dp)<=eps) exit
        end do
    end function betacf

    pure real(dp) function beta_cdf(x,a,b) result(v)
        real(dp), intent(in) :: x,a,b
        real(dp) :: bt
        if (a <= 0.0_dp .or. b <= 0.0_dp) then
            v=nan_dp(); return
        end if
        if (x <= 0.0_dp) then
            v=0.0_dp; return
        else if (x >= 1.0_dp) then
            v=1.0_dp; return
        end if
        bt=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
        if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            v=bt*betacf(a,b,x)/a
        else
            v=1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
        end if
        v=min(1.0_dp,max(0.0_dp,v))
    end function beta_cdf

end module pdqutils_special
