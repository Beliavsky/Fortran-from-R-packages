module actuar_special
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use actuar_kinds, only: dp, pi
    implicit none
    private
    public :: beta_fn, log_beta, reg_beta, inv_reg_beta
    public :: reg_gamma_p, reg_gamma_q, normal_cdf, normal_quantile
    public :: poisson_pmf, poisson_cdf, poisson_quantile
    public :: nbinom_pmf, nbinom_cdf, nbinom_quantile
    public :: binom_pmf, binom_cdf, binom_quantile
    public :: random_normal, random_gamma, random_poisson, random_binomial
    public :: random_negative_binomial, random_beta

contains

    pure real(dp) function log_beta(a, b) result(v)
        real(dp), intent(in) :: a, b
        v = log_gamma(a) + log_gamma(b) - log_gamma(a + b)
    end function log_beta

    pure real(dp) function beta_fn(a, b) result(v)
        real(dp), intent(in) :: a, b
        v = exp(log_beta(a, b))
    end function beta_fn

    pure real(dp) function beta_cf(a, b, x) result(cf)
        real(dp), intent(in) :: a, b, x
        integer, parameter :: maxit = 300
        real(dp), parameter :: eps = 3.0e-15_dp, fpmin = 1.0e-300_dp
        integer :: m, m2
        real(dp) :: aa, c, d, del, h, qab, qam, qap
        qab = a + b; qap = a + 1.0_dp; qam = a - 1.0_dp
        c = 1.0_dp
        d = 1.0_dp - qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1.0_dp/d
        h = d
        do m = 1, maxit
            m2 = 2*m
            aa = real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
            d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d; h = h*d*c
            aa = -(a+real(m,dp))*(qab+real(m,dp))*x/((a+real(m2,dp))*(qap+real(m2,dp)))
            d = 1.0_dp + aa*d; if (abs(d) < fpmin) d = fpmin
            c = 1.0_dp + aa/c; if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d; del = d*c; h = h*del
            if (abs(del - 1.0_dp) <= eps) exit
        end do
        cf = h
    end function beta_cf

    pure real(dp) function reg_beta(x, a, b) result(v)
        real(dp), intent(in) :: x, a, b
        real(dp) :: bt
        if (a <= 0.0_dp .or. b <= 0.0_dp) then
            v = ieee_value(0.0_dp, ieee_quiet_nan)
            return
        end if
        if (x <= 0.0_dp) then
            v = 0.0_dp; return
        else if (x >= 1.0_dp) then
            v = 1.0_dp; return
        end if
        bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
        if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
            v = bt*beta_cf(a,b,x)/a
        else
            v = 1.0_dp - bt*beta_cf(b,a,1.0_dp-x)/b
        end if
        v = max(0.0_dp, min(1.0_dp, v))
    end function reg_beta

    pure real(dp) function inv_reg_beta(p, a, b) result(x)
        real(dp), intent(in) :: p, a, b
        integer :: i
        real(dp) :: lo, hi, mid
        if (p <= 0.0_dp) then
            x = 0.0_dp; return
        else if (p >= 1.0_dp) then
            x = 1.0_dp; return
        end if
        lo = 0.0_dp; hi = 1.0_dp
        do i = 1, 100
            mid = 0.5_dp*(lo+hi)
            if (reg_beta(mid,a,b) < p) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = 0.5_dp*(lo+hi)
    end function inv_reg_beta

    pure real(dp) function reg_gamma_p(a, x) result(p)
        real(dp), intent(in) :: a, x
        integer, parameter :: maxit = 400
        real(dp), parameter :: eps = 3.0e-15_dp, fpmin = 1.0e-300_dp
        integer :: n
        real(dp) :: ap, del, sum, b, c, d, h, an
        if (a <= 0.0_dp .or. x < 0.0_dp) then
            p = ieee_value(0.0_dp, ieee_quiet_nan); return
        end if
        if (x == 0.0_dp) then
            p = 0.0_dp; return
        end if
        if (x < a + 1.0_dp) then
            ap = a; sum = 1.0_dp/a; del = sum
            do n = 1, maxit
                ap = ap + 1.0_dp
                del = del*x/ap
                sum = sum + del
                if (abs(del) < abs(sum)*eps) exit
            end do
            p = sum*exp(-x+a*log(x)-log_gamma(a))
        else
            b = x + 1.0_dp - a
            c = 1.0_dp/fpmin
            d = 1.0_dp/b
            h = d
            do n = 1, maxit
                an = -real(n,dp)*(real(n,dp)-a)
                b = b + 2.0_dp
                d = an*d+b; if (abs(d)<fpmin) d=fpmin
                c = b+an/c; if (abs(c)<fpmin) c=fpmin
                d = 1.0_dp/d
                del = d*c
                h = h*del
                if (abs(del-1.0_dp)<eps) exit
            end do
            p = 1.0_dp - exp(-x+a*log(x)-log_gamma(a))*h
        end if
        p = max(0.0_dp,min(1.0_dp,p))
    end function reg_gamma_p

    pure real(dp) function reg_gamma_q(a, x) result(q)
        real(dp), intent(in) :: a, x
        q = 1.0_dp - reg_gamma_p(a,x)
    end function reg_gamma_q

    pure real(dp) function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        p = 0.5_dp*erfc(-x/sqrt(2.0_dp))
    end function normal_cdf

    pure real(dp) function normal_quantile(p) result(x)
        real(dp), intent(in) :: p
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
        real(dp) :: q, r
        if (p <= 0.0_dp) then
            x = -huge(1.0_dp); return
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp); return
        end if
        if (p < plow) then
            q = sqrt(-2.0_dp*log(p))
            x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        else if (p <= phigh) then
            q = p - 0.5_dp; r = q*q
            x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
        else
            q = sqrt(-2.0_dp*log(1.0_dp-p))
            x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
        end if
    end function normal_quantile

    pure real(dp) function poisson_pmf(k, lambda) result(p)
        integer, intent(in) :: k
        real(dp), intent(in) :: lambda
        if (k < 0 .or. lambda < 0.0_dp) then
            p = 0.0_dp
        else if (lambda == 0.0_dp) then
            p = merge(1.0_dp,0.0_dp,k==0)
        else
            p = exp(-lambda + real(k,dp)*log(lambda)-log_gamma(real(k+1,dp)))
        end if
    end function poisson_pmf

    pure real(dp) function poisson_cdf(k, lambda) result(p)
        integer, intent(in) :: k
        real(dp), intent(in) :: lambda
        integer :: j
        real(dp) :: term
        if (k < 0) then; p=0.0_dp; return; end if
        if (lambda == 0.0_dp) then; p=1.0_dp; return; end if
        term = exp(-lambda); p = term
        do j=1,k
            term = term*lambda/real(j,dp); p=p+term
        end do
        p=min(1.0_dp,p)
    end function poisson_cdf

    pure integer function poisson_quantile(prob, lambda) result(k)
        real(dp), intent(in) :: prob, lambda
        real(dp) :: c, term
        if (prob<=0.0_dp) then; k=0; return; end if
        term=exp(-lambda); c=term; k=0
        do while (c < prob .and. k < 100000)
            k=k+1; term=term*lambda/real(k,dp); c=c+term
        end do
    end function poisson_quantile

    pure real(dp) function nbinom_pmf(k, size, prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob
        if(k<0 .or. size<=0.0_dp .or. prob<=0.0_dp .or. prob>1.0_dp) then
            p=0.0_dp
        else
            p=exp(log_gamma(real(k,dp)+size)-log_gamma(size)-log_gamma(real(k+1,dp)) &
                + size*log(prob)+real(k,dp)*log(1.0_dp-prob))
        end if
    end function nbinom_pmf

    pure real(dp) function nbinom_cdf(k,size,prob) result(p)
        integer,intent(in)::k
        real(dp),intent(in)::size,prob
        if(k<0) then; p=0.0_dp; else; p=reg_beta(prob,size,real(k+1,dp)); end if
    end function nbinom_cdf

    pure integer function nbinom_quantile(q,size,prob) result(k)
        real(dp),intent(in)::q,size,prob
        if(q<=0.0_dp) then;k=0;return;end if
        k=0
        do while(nbinom_cdf(k,size,prob)<q .and. k<100000); k=k+1; end do
    end function nbinom_quantile

    pure real(dp) function binom_pmf(k,n,prob) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob
        if(k<0 .or. k>n .or. prob<0.0_dp .or. prob>1.0_dp) then;p=0.0_dp;return;end if
        if(prob==0.0_dp) then;p=merge(1.0_dp,0.0_dp,k==0);return;end if
        if(prob==1.0_dp) then;p=merge(1.0_dp,0.0_dp,k==n);return;end if
        p=exp(log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp)) &
              +real(k,dp)*log(prob)+real(n-k,dp)*log(1.0_dp-prob))
    end function binom_pmf

    pure real(dp) function binom_cdf(k,n,prob) result(p)
        integer,intent(in)::k,n
        real(dp),intent(in)::prob
        integer::j
        if(k<0) then;p=0.0_dp;return;end if
        if(k>=n) then;p=1.0_dp;return;end if
        p=0.0_dp
        do j=0,k; p=p+binom_pmf(j,n,prob); end do
    end function binom_cdf

    pure integer function binom_quantile(q,n,prob) result(k)
        real(dp),intent(in)::q,prob
        integer,intent(in)::n
        k=0; do while(k<n .and. binom_cdf(k,n,prob)<q); k=k+1; end do
    end function binom_quantile

    real(dp) function random_normal() result(z)
        real(dp) :: u1,u2
        call random_number(u1); call random_number(u2)
        u1=max(u1,tiny(1.0_dp)); z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
    end function random_normal

    recursive real(dp) function random_gamma(shape, scale) result(x)
        real(dp),intent(in)::shape,scale
        real(dp)::d,c,z,u,v
        if(shape<=0.0_dp .or. scale<=0.0_dp) then;x=ieee_value(0.0_dp, ieee_quiet_nan);return;end if
        if(shape<1.0_dp) then
            call random_number(u); x=random_gamma(shape+1.0_dp,scale)*u**(1.0_dp/shape); return
        end if
        d=shape-1.0_dp/3.0_dp; c=1.0_dp/sqrt(9.0_dp*d)
        do
            z=random_normal(); v=(1.0_dp+c*z)**3
            if(v<=0.0_dp) cycle
            call random_number(u)
            if(u<1.0_dp-0.0331_dp*z**4) exit
            if(log(u)<0.5_dp*z*z+d*(1.0_dp-v+log(v))) exit
        end do
        x=scale*d*v
    end function random_gamma

    integer function random_poisson(lambda) result(k)
        real(dp),intent(in)::lambda
        real(dp)::l,p,u
        if(lambda<30.0_dp) then
            l=exp(-lambda); p=1.0_dp; k=-1
            do; k=k+1; call random_number(u); p=p*u; if(p<=l) exit; end do
        else
            k=max(0,nint(lambda+sqrt(lambda)*random_normal()))
        end if
    end function random_poisson

    integer function random_binomial(n,prob) result(k)
        integer,intent(in)::n
        real(dp),intent(in)::prob
        integer::i
        real(dp)::u
        k=0; do i=1,n; call random_number(u); if(u<prob) k=k+1; end do
    end function random_binomial

    integer function random_negative_binomial(size,prob) result(k)
        real(dp),intent(in)::size,prob
        real(dp)::lambda
        lambda=random_gamma(size,(1.0_dp-prob)/prob)
        k=random_poisson(lambda)
    end function random_negative_binomial

    real(dp) function random_beta(a,b) result(x)
        real(dp),intent(in)::a,b
        real(dp)::g1,g2
        g1=random_gamma(a,1.0_dp); g2=random_gamma(b,1.0_dp); x=g1/(g1+g2)
    end function random_beta

end module actuar_special
