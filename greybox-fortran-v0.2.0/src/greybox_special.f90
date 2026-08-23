module greybox_special
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use greybox_kinds, only: dp, pi
    implicit none
    private
    public :: normal_pdf, normal_cdf, normal_quantile, normal_rng
    public :: gamma_p, gamma_quantile, gamma_rng, beta_inc
    public :: student_t_pdf, student_t_cdf, student_t_quantile
    public :: logistic_pdf, logistic_cdf, logistic_quantile
    public :: log1pexp, safe_log, nan_dp, inf_dp

contains

    pure elemental real(dp) function nan_dp() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_dp

    pure elemental real(dp) function inf_dp() result(x)
        x = ieee_value(0.0_dp, ieee_positive_inf)
    end function inf_dp

    pure elemental real(dp) function safe_log(x) result(y)
        real(dp), intent(in) :: x
        if (x > 0.0_dp) then
            y = log(x)
        else if (x <= 0.0_dp) then
            y = -inf_dp()
        else
            y = nan_dp()
        end if
    end function safe_log

    pure elemental real(dp) function log1pexp(x) result(y)
        real(dp), intent(in) :: x
        if (x > 0.0_dp) then
            y = x + log(1.0_dp + exp(-x))
        else
            y = log(1.0_dp + exp(x))
        end if
    end function log1pexp

    pure elemental real(dp) function normal_pdf(x, mu, sigma) result(f)
        real(dp), intent(in) :: x, mu, sigma
        real(dp) :: z
        if (sigma <= 0.0_dp) then
            f = nan_dp()
            return
        end if
        z = (x - mu) / sigma
        f = exp(-0.5_dp*z*z) / (sqrt(2.0_dp*pi)*sigma)
    end function normal_pdf

    pure elemental real(dp) function normal_cdf(x, mu, sigma) result(p)
        real(dp), intent(in) :: x, mu, sigma
        if (sigma <= 0.0_dp) then
            p = nan_dp()
        else
            p = 0.5_dp * erfc(-(x-mu)/(sigma*sqrt(2.0_dp)))
        end if
    end function normal_cdf

    pure elemental real(dp) function normal_quantile(p, mu, sigma) result(x)
        real(dp), intent(in) :: p, mu, sigma
        real(dp), parameter :: a1=-3.969683028665376d1, a2=2.209460984245205d2
        real(dp), parameter :: a3=-2.759285104469687d2, a4=1.383577518672690d2
        real(dp), parameter :: a5=-3.066479806614716d1, a6=2.506628277459239d0
        real(dp), parameter :: b1=-5.447609879822406d1, b2=1.615858368580409d2
        real(dp), parameter :: b3=-1.556989798598866d2, b4=6.680131188771972d1
        real(dp), parameter :: b5=-1.328068155288572d1
        real(dp), parameter :: c1=-7.784894002430293d-3, c2=-3.223964580411365d-1
        real(dp), parameter :: c3=-2.400758277161838d0, c4=-2.549732539343734d0
        real(dp), parameter :: c5=4.374664141464968d0, c6=2.938163982698783d0
        real(dp), parameter :: d1=7.784695709041462d-3, d2=3.224671290700398d-1
        real(dp), parameter :: d3=2.445134137142996d0, d4=3.754408661907416d0
        real(dp), parameter :: plow=0.02425_dp, phigh=1.0_dp-plow
        real(dp) :: q, r, z
        if (sigma <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = nan_dp(); return
        else if (p <= 0.0_dp) then
            x = -inf_dp(); return
        else if (p >= 1.0_dp) then
            x = inf_dp(); return
        end if
        if (p < plow) then
            q = sqrt(-2.0_dp*log(p))
            z = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        else if (p <= phigh) then
            q = p - 0.5_dp
            r = q*q
            z = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
                (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
        else
            q = sqrt(-2.0_dp*log(1.0_dp-p))
            z = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
                 ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        end if
        ! One Halley refinement.
        q = normal_cdf(z,0.0_dp,1.0_dp) - p
        z = z - q / max(normal_pdf(z,0.0_dp,1.0_dp), tiny(1.0_dp))
        x = mu + sigma*z
    end function normal_quantile

    real(dp) function normal_rng(mu, sigma) result(x)
        real(dp), intent(in) :: mu, sigma
        real(dp) :: u1, u2
        if (sigma < 0.0_dp) then
            x = nan_dp(); return
        end if
        call random_number(u1); call random_number(u2)
        u1 = max(u1, tiny(1.0_dp))
        x = mu + sigma*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
    end function normal_rng

    pure elemental real(dp) function gamma_p(a, x) result(p)
        real(dp), intent(in) :: a, x
        integer, parameter :: itmax=500
        real(dp), parameter :: eps=2.0e-14_dp, fpmin=1.0e-300_dp
        integer :: n
        real(dp) :: ap, del, sum, b, c, d, h, an
        if (a <= 0.0_dp .or. x < 0.0_dp) then
            p = nan_dp(); return
        else if (x <= 0.0_dp) then
            p = 0.0_dp; return
        end if
        if (x < a + 1.0_dp) then
            ap = a; sum = 1.0_dp/a; del = sum
            do n=1,itmax
                ap = ap + 1.0_dp
                del = del*x/ap
                sum = sum + del
                if (abs(del) <= abs(sum)*eps) exit
            end do
            p = sum*exp(-x+a*log(x)-log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp/fpmin
            d = 1.0_dp/b
            h = d
            do n=1,itmax
                an = -real(n,dp)*(real(n,dp)-a)
                b = b + 2.0_dp
                d = an*d + b
                if (abs(d) < fpmin) d = fpmin
                c = b + an/c
                if (abs(c) < fpmin) c = fpmin
                d = 1.0_dp/d
                del = d*c
                h = h*del
                if (abs(del-1.0_dp) <= eps) exit
            end do
            p = 1.0_dp - exp(-x+a*log(x)-log_gamma(a))*h
        end if
        p = min(1.0_dp,max(0.0_dp,p))
    end function gamma_p

    pure elemental real(dp) function gamma_quantile(p, a, scale) result(x)
        real(dp), intent(in) :: p, a, scale
        real(dp) :: lo, hi, mid, z, f, deriv
        integer :: i
        if (a <= 0.0_dp .or. scale <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = nan_dp(); return
        else if (p <= 0.0_dp) then
            x = 0.0_dp; return
        else if (p >= 1.0_dp) then
            x = inf_dp(); return
        end if
        z = normal_quantile(p,0.0_dp,1.0_dp)
        if (a > 0.2_dp) then
            x = a*max(1.0e-6_dp, 1.0_dp - 1.0_dp/(9.0_dp*a) + z/(3.0_dp*sqrt(a)))**3
        else
            x = max(1.0e-12_dp, (p*exp(log_gamma(a+1.0_dp)))**(1.0_dp/a))
        end if
        do i=1,8
            f = gamma_p(a,x)-p
            deriv = exp((a-1.0_dp)*log(max(x,tiny(1.0_dp)))-x-log_gamma(a))
            if (deriv <= tiny(1.0_dp)) exit
            mid = x - f/deriv
            if (mid <= 0.0_dp .or. .not.(mid < huge(1.0_dp))) exit
            if (abs(mid-x) <= 1.0e-12_dp*max(1.0_dp,x)) then
                x = mid; exit
            end if
            x = mid
        end do
        lo = 0.0_dp; hi = max(1.0_dp,x)
        do while (gamma_p(a,hi) < p)
            hi = 2.0_dp*hi
            if (hi > huge(1.0_dp)/4.0_dp) exit
        end do
        do i=1,100
            mid = 0.5_dp*(lo+hi)
            if (gamma_p(a,mid) < p) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = 0.5_dp*(lo+hi)*scale
    end function gamma_quantile

    recursive real(dp) function gamma_rng(shape, scale) result(x)
        real(dp), intent(in) :: shape, scale
        real(dp) :: d, c, z, u, v
        if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
            x = nan_dp(); return
        end if
        if (shape < 1.0_dp) then
            call random_number(u)
            x = gamma_rng(shape+1.0_dp,scale)*u**(1.0_dp/shape)
            return
        end if
        d = shape - 1.0_dp/3.0_dp
        c = 1.0_dp/sqrt(9.0_dp*d)
        do
            z = normal_rng(0.0_dp,1.0_dp)
            v = (1.0_dp+c*z)**3
            if (v <= 0.0_dp) cycle
            call random_number(u)
            if (u < 1.0_dp-0.0331_dp*z**4) exit
            if (log(u) < 0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
        end do
        x = scale*d*v
    end function gamma_rng

    pure elemental real(dp) function beta_inc(a,b,x) result(v)
        real(dp), intent(in) :: a,b,x
        real(dp) :: bt
        if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
            v = nan_dp(); return
        else if (x <= 0.0_dp) then
            v = 0.0_dp; return
        else if (x >= 1.0_dp) then
            v = 1.0_dp; return
        end if
        bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
        if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            v = bt*betacf(a,b,x)/a
        else
            v = 1.0_dp-bt*betacf(b,a,1.0_dp-x)/b
        end if
        v = min(1.0_dp,max(0.0_dp,v))
    end function beta_inc

    pure elemental real(dp) function betacf(a,b,x) result(h)
        real(dp), intent(in) :: a,b,x
        integer, parameter :: maxit=300
        real(dp), parameter :: eps=3.0e-14_dp, fpmin=1.0e-300_dp
        integer :: m, m2
        real(dp) :: aa,c,d,del,qab,qam,qap
        qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
        c=1.0_dp; d=1.0_dp-qab*x/qap
        if (abs(d)<fpmin) d=fpmin
        d=1.0_dp/d; h=d
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
            if(abs(del-1.0_dp)<eps)exit
        end do
    end function betacf

    pure elemental real(dp) function student_t_pdf(x, nu) result(f)
        real(dp), intent(in) :: x, nu
        if (nu <= 0.0_dp) then
            f = nan_dp()
        else
            f = exp(log_gamma((nu+1.0_dp)/2.0_dp)-log_gamma(nu/2.0_dp)) / &
                sqrt(nu*pi) * (1.0_dp+x*x/nu)**(-(nu+1.0_dp)/2.0_dp)
        end if
    end function student_t_pdf

    pure elemental real(dp) function student_t_cdf(x, nu) result(p)
        real(dp), intent(in) :: x, nu
        real(dp) :: ib
        if (nu <= 0.0_dp) then
            p = nan_dp(); return
        else if (x <= 0.0_dp) then
            p = 0.5_dp; return
        end if
        ib = beta_inc(nu/2.0_dp,0.5_dp,nu/(nu+x*x))
        if (x > 0.0_dp) then
            p = 1.0_dp-0.5_dp*ib
        else
            p = 0.5_dp*ib
        end if
    end function student_t_cdf

    pure elemental real(dp) function student_t_quantile(p, nu) result(x)
        real(dp), intent(in) :: p, nu
        real(dp) :: lo,hi,mid
        integer :: i
        if (nu <= 0.0_dp .or. p < 0.0_dp .or. p > 1.0_dp) then
            x = nan_dp(); return
        else if (p<=0.0_dp) then
            x=-inf_dp(); return
        else if (p>=1.0_dp) then
            x=inf_dp(); return
        else if (abs(p-0.5_dp)<=epsilon(1.0_dp)) then
            x=0.0_dp; return
        end if
        lo=-1.0_dp; hi=1.0_dp
        do while(student_t_cdf(lo,nu)>p); lo=2.0_dp*lo; end do
        do while(student_t_cdf(hi,nu)<p); hi=2.0_dp*hi; end do
        do i=1,120
            mid=0.5_dp*(lo+hi)
            if(student_t_cdf(mid,nu)<p)then; lo=mid; else; hi=mid; end if
        end do
        x=0.5_dp*(lo+hi)
    end function student_t_quantile

    pure elemental real(dp) function logistic_pdf(x, mu, scale) result(f)
        real(dp), intent(in) :: x,mu,scale
        real(dp) :: z
        if(scale<=0.0_dp)then; f=nan_dp(); return; end if
        z=(x-mu)/scale
        if(z>=0.0_dp)then
            f=exp(-z)/(scale*(1.0_dp+exp(-z))**2)
        else
            f=exp(z)/(scale*(1.0_dp+exp(z))**2)
        end if
    end function logistic_pdf

    pure elemental real(dp) function logistic_cdf(x,mu,scale) result(p)
        real(dp),intent(in)::x,mu,scale
        real(dp)::z
        if(scale<=0.0_dp)then;p=nan_dp();return;end if
        z=(x-mu)/scale
        if(z>=0.0_dp)then;p=1.0_dp/(1.0_dp+exp(-z));else;p=exp(z)/(1.0_dp+exp(z));end if
    end function logistic_cdf

    pure elemental real(dp) function logistic_quantile(p,mu,scale) result(x)
        real(dp),intent(in)::p,mu,scale
        if(scale<=0.0_dp.or.p<0.0_dp.or.p>1.0_dp)then;x=nan_dp();return;end if
        if (p <= 0.0_dp) then
            x = -inf_dp()
        else if (p >= 1.0_dp) then
            x = inf_dp()
        else
            x = mu + scale*log(p/(1.0_dp-p))
        end if
    end function logistic_quantile

end module greybox_special
