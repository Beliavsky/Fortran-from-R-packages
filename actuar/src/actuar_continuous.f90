module actuar_continuous
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use actuar_kinds, only: dp, pi
    use actuar_special, only: beta_fn, log_beta, reg_beta, inv_reg_beta, &
                              reg_gamma_p, reg_gamma_q, normal_cdf, normal_quantile, &
                              random_beta, random_gamma, random_normal
    implicit none
    private
    public :: dpareto, ppareto, qpareto, rpareto, mpareto, levpareto
    public :: dpareto1, ppareto1, qpareto1, rpareto1, mpareto1, levpareto1
    public :: dpareto2, ppareto2, qpareto2, rpareto2
    public :: dpareto3, ppareto3, qpareto3, rpareto3
    public :: dpareto4, ppareto4, qpareto4, rpareto4
    public :: dburr, pburr, qburr, rburr, mburr
    public :: dinvburr, pinvburr, qinvburr, rinvburr, minvburr
    public :: dllogis, pllogis, qllogis, rllogis, mllogis
    public :: dparalogis, pparalogis, qparalogis, rparalogis
    public :: dinvparalogis, pinvparalogis, qinvparalogis, rinvparalogis
    public :: dgenpareto, pgenpareto, qgenpareto, rgenpareto, mgenpareto
    public :: dgenbeta, pgenbeta, qgenbeta, rgenbeta, mgenbeta
    public :: dinvexp, pinvexp, qinvexp, rinvexp, minvexp
    public :: dinvgamma, pinvgamma, qinvgamma, rinvgamma, minvgamma
    public :: dinvweibull, pinvweibull, qinvweibull, rinvweibull, minvweibull
    public :: dgumbel, pgumbel, qgumbel, rgumbel, mgumbel, mgfgumbel
    public :: dlgamma, plgamma, qlgamma, rlgamma, mlgamma
    public :: dtrgamma, ptrgamma, qtrgamma, rtrgamma, mtrgamma
    public :: dtrbeta, ptrbeta, qtrbeta, rtrbeta, mtrbeta

contains

    pure real(dp) function dpareto(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale
        if(x<0.0_dp .or. shape<=0.0_dp .or. scale<=0.0_dp) then; f=0.0_dp; return; end if
        f=shape/scale*(1.0_dp+x/scale)**(-shape-1.0_dp)
    end function dpareto
    pure real(dp) function ppareto(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=1.0_dp-(1.0_dp+x/scale)**(-shape);end if
    end function ppareto
    pure real(dp) function qpareto(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale
        if(p<=0.0_dp) then;x=0.0_dp;else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=scale*((1.0_dp-p)**(-1.0_dp/shape)-1.0_dp);end if
    end function qpareto
    real(dp) function rpareto(shape,scale) result(x)
        real(dp),intent(in)::shape,scale; real(dp)::u
        call random_number(u); x=qpareto(u,shape,scale)
    end function rpareto
    pure real(dp) function mpareto(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        if(order>=shape .or. order<=-1.0_dp) then;m=huge(1.0_dp);else
            m=scale**order*gamma(1.0_dp+order)*gamma(shape-order)/gamma(shape); end if
    end function mpareto
    pure real(dp) function levpareto(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::u,a,b
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        if(order<=-1.0_dp) then;m=huge(1.0_dp);return;end if
        u=limit/(limit+scale); a=1.0_dp+order; b=shape-order
        if(b>0.0_dp) then
            m=scale**order*beta_fn(a,b)*reg_beta(u,a,b)*shape + limit**order*(1.0_dp-ppareto(limit,shape,scale))
        else
            m=numeric_limited_moment_pareto(limit,shape,scale,order)
        end if
    end function levpareto
    pure real(dp) function numeric_limited_moment_pareto(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        integer,parameter::n=1000
        integer::i; real(dp)::h,x,s
        h=limit/real(n,dp); s=0.0_dp
        do i=0,n
            x=h*real(i,dp)
            if(i==0 .or. i==n) then;s=s+0.5_dp*x**order*dpareto(x,shape,scale);else;s=s+x**order*dpareto(x,shape,scale);end if
        end do
        m=h*s+limit**order*(1.0_dp-ppareto(limit,shape,scale))
    end function numeric_limited_moment_pareto

    pure real(dp) function dpareto1(x,shape,xmin) result(f)
        real(dp),intent(in)::x,shape,xmin
        if(x<xmin .or. shape<=0.0_dp .or. xmin<=0.0_dp) then;f=0.0_dp;else;f=shape*xmin**shape/x**(shape+1.0_dp);end if
    end function dpareto1
    pure real(dp) function ppareto1(x,shape,xmin) result(p)
        real(dp),intent(in)::x,shape,xmin
        if(x<xmin) then;p=0.0_dp;else;p=1.0_dp-(xmin/x)**shape;end if
    end function ppareto1
    pure real(dp) function qpareto1(p,shape,xmin) result(x)
        real(dp),intent(in)::p,shape,xmin
        if(p<=0.0_dp) then;x=xmin;else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=xmin/(1.0_dp-p)**(1.0_dp/shape);end if
    end function qpareto1
    real(dp) function rpareto1(shape,xmin) result(x)
        real(dp),intent(in)::shape,xmin; real(dp)::u
        call random_number(u);x=qpareto1(u,shape,xmin)
    end function rpareto1
    pure real(dp) function mpareto1(order,shape,xmin) result(m)
        real(dp),intent(in)::order,shape,xmin
        if(order>=shape) then;m=huge(1.0_dp);else;m=shape*xmin**order/(shape-order);end if
    end function mpareto1
    pure real(dp) function levpareto1(limit,shape,xmin,order) result(m)
        real(dp),intent(in)::limit,shape,xmin,order
        if(limit<=xmin) then;m=limit**order;return;end if
        if(abs(shape-order)<1.0e-12_dp) then
            m=shape*xmin**shape*log(limit/xmin)+limit**order*(xmin/limit)**shape
        else
            m=shape*xmin**shape*(limit**(order-shape)-xmin**(order-shape))/(order-shape) &
              +limit**order*(xmin/limit)**shape
        end if
    end function levpareto1

    pure real(dp) function dpareto2(x,xmin,shape,scale) result(f)
        real(dp),intent(in)::x,xmin,shape,scale
        f=merge(dpareto(x-xmin,shape,scale),0.0_dp,x>=xmin)
    end function dpareto2
    pure real(dp) function ppareto2(x,xmin,shape,scale) result(p)
        real(dp),intent(in)::x,xmin,shape,scale
        if(x<xmin) then;p=0.0_dp;else;p=ppareto(x-xmin,shape,scale);end if
    end function ppareto2
    pure real(dp) function qpareto2(p,xmin,shape,scale) result(x)
        real(dp),intent(in)::p,xmin,shape,scale;x=xmin+qpareto(p,shape,scale)
    end function qpareto2
    real(dp) function rpareto2(xmin,shape,scale) result(x)
        real(dp),intent(in)::xmin,shape,scale;x=xmin+rpareto(shape,scale)
    end function rpareto2

    pure real(dp) function dllogis(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale; real(dp)::v
        if(x<=0.0_dp .or. shape<=0.0_dp .or. scale<=0.0_dp) then;f=0.0_dp;return;end if
        v=(x/scale)**shape;f=shape*v/(x*(1.0_dp+v)**2)
    end function dllogis
    pure real(dp) function pllogis(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale; real(dp)::v
        if(x<=0.0_dp) then;p=0.0_dp;else;v=(x/scale)**shape;p=v/(1.0_dp+v);end if
    end function pllogis
    pure real(dp) function qllogis(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale
        if(p<=0.0_dp) then;x=0.0_dp;else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=scale*(p/(1.0_dp-p))**(1.0_dp/shape);end if
    end function qllogis
    real(dp) function rllogis(shape,scale) result(x)
        real(dp),intent(in)::shape,scale;real(dp)::u;call random_number(u);x=qllogis(u,shape,scale)
    end function rllogis
    pure real(dp) function mllogis(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale; real(dp)::t
        if(order<=-shape .or. order>=shape) then
            m=huge(1.0_dp)
        else
            t=order/shape
            m=scale**order*gamma(1.0_dp+t)*gamma(1.0_dp-t)
        end if
    end function mllogis

    pure real(dp) function dpareto3(x,xmin,shape,scale) result(f)
        real(dp),intent(in)::x,xmin,shape,scale
        if(x<xmin) then;f=0.0_dp;else;f=dllogis(x-xmin,shape,scale);end if
    end function dpareto3
    pure real(dp) function ppareto3(x,xmin,shape,scale) result(p)
        real(dp),intent(in)::x,xmin,shape,scale
        if(x<xmin) then;p=0.0_dp;else;p=pllogis(x-xmin,shape,scale);end if
    end function ppareto3
    pure real(dp) function qpareto3(p,xmin,shape,scale) result(x)
        real(dp),intent(in)::p,xmin,shape,scale;x=xmin+qllogis(p,shape,scale)
    end function qpareto3
    real(dp) function rpareto3(xmin,shape,scale) result(x)
        real(dp),intent(in)::xmin,shape,scale;x=xmin+rllogis(shape,scale)
    end function rpareto3

    pure real(dp) function dburr(x,shape1,shape2,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,scale; real(dp)::v
        if(x<=0.0_dp .or. min(shape1,shape2,scale)<=0.0_dp) then;f=0.0_dp;return;end if
        v=(x/scale)**shape2
        f=shape1*shape2/scale*(x/scale)**(shape2-1.0_dp)*(1.0_dp+v)**(-shape1-1.0_dp)
    end function dburr
    pure real(dp) function pburr(x,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=1.0_dp-(1.0_dp+(x/scale)**shape2)**(-shape1);end if
    end function pburr
    pure real(dp) function qburr(p,shape1,shape2,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,scale
        if(p<=0.0_dp) then
            x=0.0_dp
        else if(p>=1.0_dp) then
            x=huge(1.0_dp)
        else
            x=scale*((1.0_dp-p)**(-1.0_dp/shape1)-1.0_dp)**(1.0_dp/shape2)
        end if
    end function qburr
    real(dp) function rburr(shape1,shape2,scale) result(x)
        real(dp),intent(in)::shape1,shape2,scale;real(dp)::u;call random_number(u);x=qburr(u,shape1,shape2,scale)
    end function rburr
    pure real(dp) function mburr(order,shape1,shape2,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,scale; real(dp)::t
        t=order/shape2
        if(t<=-1.0_dp .or. t>=shape1) then;m=huge(1.0_dp);else;m=scale**order*gamma(1.0_dp+t)*gamma(shape1-t)/gamma(shape1);end if
    end function mburr

    pure real(dp) function dpareto4(x,xmin,shape1,shape2,scale) result(f)
        real(dp),intent(in)::x,xmin,shape1,shape2,scale
        if(x<xmin) then;f=0.0_dp;else;f=dburr(x-xmin,shape1,shape2,scale);end if
    end function dpareto4
    pure real(dp) function ppareto4(x,xmin,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,xmin,shape1,shape2,scale
        if(x<xmin) then;p=0.0_dp;else;p=pburr(x-xmin,shape1,shape2,scale);end if
    end function ppareto4
    pure real(dp) function qpareto4(p,xmin,shape1,shape2,scale) result(x)
        real(dp),intent(in)::p,xmin,shape1,shape2,scale;x=xmin+qburr(p,shape1,shape2,scale)
    end function qpareto4
    real(dp) function rpareto4(xmin,shape1,shape2,scale) result(x)
        real(dp),intent(in)::xmin,shape1,shape2,scale;x=xmin+rburr(shape1,shape2,scale)
    end function rpareto4

    pure real(dp) function dinvburr(x,shape1,shape2,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,scale; real(dp)::v,u
        if(x<=0.0_dp .or. min(shape1,shape2,scale)<=0.0_dp) then;f=0.0_dp;return;end if
        v=(scale/x)**shape2;u=1.0_dp/(1.0_dp+v)
        f=shape1*shape2*u**shape1*(1.0_dp-u)/x
    end function dinvburr
    pure real(dp) function pinvburr(x,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=(1.0_dp+(scale/x)**shape2)**(-shape1);end if
    end function pinvburr
    pure real(dp) function qinvburr(p,shape1,shape2,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,scale
        if(p<=0.0_dp) then
            x=0.0_dp
        else if(p>=1.0_dp) then
            x=huge(1.0_dp)
        else
            x=scale*(p**(-1.0_dp/shape1)-1.0_dp)**(-1.0_dp/shape2)
        end if
    end function qinvburr
    real(dp) function rinvburr(shape1,shape2,scale) result(x)
        real(dp),intent(in)::shape1,shape2,scale;real(dp)::u;call random_number(u);x=qinvburr(u,shape1,shape2,scale)
    end function rinvburr
    pure real(dp) function minvburr(order,shape1,shape2,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,scale; real(dp)::t
        t=order/shape2
        if(t>=1.0_dp .or. shape1+t<=0.0_dp) then
            m=huge(1.0_dp)
        else
            m=scale**order*gamma(shape1+t)*gamma(1.0_dp-t)/gamma(shape1)
        end if
    end function minvburr

    pure real(dp) function dparalogis(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale;f=dburr(x,shape,shape,scale)
    end function dparalogis
    pure real(dp) function pparalogis(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale;p=pburr(x,shape,shape,scale)
    end function pparalogis
    pure real(dp) function qparalogis(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale;x=qburr(p,shape,shape,scale)
    end function qparalogis
    real(dp) function rparalogis(shape,scale) result(x)
        real(dp),intent(in)::shape,scale;x=rburr(shape,shape,scale)
    end function rparalogis
    pure real(dp) function dinvparalogis(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale;f=dinvburr(x,shape,shape,scale)
    end function dinvparalogis
    pure real(dp) function pinvparalogis(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale;p=pinvburr(x,shape,shape,scale)
    end function pinvparalogis
    pure real(dp) function qinvparalogis(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale;x=qinvburr(p,shape,shape,scale)
    end function qinvparalogis
    real(dp) function rinvparalogis(shape,scale) result(x)
        real(dp),intent(in)::shape,scale;x=rinvburr(shape,shape,scale)
    end function rinvparalogis

    pure real(dp) function dgenpareto(x,shape1,shape2,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,scale; real(dp)::u
        if(x<=0.0_dp .or. min(shape1,shape2,scale)<=0.0_dp) then;f=0.0_dp;return;end if
        u=x/(x+scale)
        f=u**shape2*(1.0_dp-u)**shape1/(x*beta_fn(shape2,shape1))
    end function dgenpareto
    pure real(dp) function pgenpareto(x,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=reg_beta(x/(x+scale),shape2,shape1);end if
    end function pgenpareto
    pure real(dp) function qgenpareto(p,shape1,shape2,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,scale;real(dp)::u
        u=inv_reg_beta(p,shape2,shape1);if(u>=1.0_dp) then;x=huge(1.0_dp);else;x=scale*u/(1.0_dp-u);end if
    end function qgenpareto
    real(dp) function rgenpareto(shape1,shape2,scale) result(x)
        real(dp),intent(in)::shape1,shape2,scale;real(dp)::u;u=random_beta(shape2,shape1);x=scale*u/(1.0_dp-u)
    end function rgenpareto
    pure real(dp) function mgenpareto(order,shape1,shape2,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,scale
        if(order<=-shape2 .or. order>=shape1) then
            m=huge(1.0_dp)
        else
            m=scale**order*beta_fn(shape1-order,shape2+order)/beta_fn(shape1,shape2)
        end if
    end function mgenpareto

    pure real(dp) function dgenbeta(x,shape1,shape2,shape3,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,shape3,scale;real(dp)::u
        if(x<=0.0_dp .or. x>=scale .or. min(shape1,shape2,shape3,scale)<=0.0_dp) then;f=0.0_dp;return;end if
        u=(x/scale)**shape3
        f=shape3*u**shape1*(1.0_dp-u)**(shape2-1.0_dp)/(x*beta_fn(shape1,shape2))
    end function dgenbeta
    pure real(dp) function pgenbeta(x,shape1,shape2,shape3,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,shape3,scale
        if(x<=0.0_dp) then;p=0.0_dp;else if(x>=scale) then;p=1.0_dp;else;p=reg_beta((x/scale)**shape3,shape1,shape2);end if
    end function pgenbeta
    pure real(dp) function qgenbeta(p,shape1,shape2,shape3,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,shape3,scale;x=scale*inv_reg_beta(p,shape1,shape2)**(1.0_dp/shape3)
    end function qgenbeta
    real(dp) function rgenbeta(shape1,shape2,shape3,scale) result(x)
        real(dp),intent(in)::shape1,shape2,shape3,scale;x=scale*random_beta(shape1,shape2)**(1.0_dp/shape3)
    end function rgenbeta
    pure real(dp) function mgenbeta(order,shape1,shape2,shape3,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,shape3,scale;real(dp)::t
        t=order/shape3
        if(shape1+t<=0.0_dp) then;m=huge(1.0_dp);else;m=scale**order*beta_fn(shape1+t,shape2)/beta_fn(shape1,shape2);end if
    end function mgenbeta

    pure real(dp) function dinvexp(x,scale) result(f)
        real(dp),intent(in)::x,scale
        if(x<=0.0_dp .or. scale<=0.0_dp) then;f=0.0_dp;else;f=scale*exp(-scale/x)/(x*x);end if
    end function dinvexp
    pure real(dp) function pinvexp(x,scale) result(p)
        real(dp),intent(in)::x,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=exp(-scale/x);end if
    end function pinvexp
    pure real(dp) function qinvexp(p,scale) result(x)
        real(dp),intent(in)::p,scale
        if(p<=0.0_dp) then;x=0.0_dp;else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=-scale/log(p);end if
    end function qinvexp
    real(dp) function rinvexp(scale) result(x)
        real(dp),intent(in)::scale;real(dp)::u;call random_number(u);x=qinvexp(u,scale)
    end function rinvexp
    pure real(dp) function minvexp(order,scale) result(m)
        real(dp),intent(in)::order,scale
        if(order>=1.0_dp) then;m=huge(1.0_dp);else;m=scale**order*gamma(1.0_dp-order);end if
    end function minvexp

    pure real(dp) function dinvgamma(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale
        if(x<=0.0_dp .or. min(shape,scale)<=0.0_dp) then
            f=0.0_dp
        else
            f=exp(shape*log(scale)-log_gamma(shape)-(shape+1.0_dp)*log(x)-scale/x)
        end if
    end function dinvgamma
    pure real(dp) function pinvgamma(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=reg_gamma_q(shape,scale/x);end if
    end function pinvgamma
    pure real(dp) function qinvgamma(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale;integer::i;real(dp)::lo,hi,mid
        if(p<=0.0_dp) then;x=0.0_dp;return;end if
        if(p>=1.0_dp) then;x=huge(1.0_dp);return;end if
        lo=tiny(1.0_dp);hi=max(scale,1.0_dp)
        do while(pinvgamma(hi,shape,scale)<p);hi=2.0_dp*hi;end do
        do i=1,100;mid=0.5_dp*(lo+hi);if(pinvgamma(mid,shape,scale)<p) then;lo=mid;else;hi=mid;end if;end do
        x=0.5_dp*(lo+hi)
    end function qinvgamma
    real(dp) function rinvgamma(shape,scale) result(x)
        real(dp),intent(in)::shape,scale;x=1.0_dp/random_gamma(shape,1.0_dp/scale)
    end function rinvgamma
    pure real(dp) function minvgamma(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        if(order>=shape) then;m=huge(1.0_dp);else;m=scale**order*gamma(shape-order)/gamma(shape);end if
    end function minvgamma

    pure real(dp) function dinvweibull(x,shape,scale) result(f)
        real(dp),intent(in)::x,shape,scale;real(dp)::z
        if(x<=0.0_dp .or. min(shape,scale)<=0.0_dp) then
            f=0.0_dp
        else
            z=(scale/x)**shape
            f=shape*scale**shape*x**(-shape-1.0_dp)*exp(-z)
        end if
    end function dinvweibull
    pure real(dp) function pinvweibull(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=exp(-(scale/x)**shape);end if
    end function pinvweibull
    pure real(dp) function qinvweibull(p,shape,scale) result(x)
        real(dp),intent(in)::p,shape,scale
        if(p<=0.0_dp) then;x=0.0_dp;else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=scale/(-log(p))**(1.0_dp/shape);end if
    end function qinvweibull
    real(dp) function rinvweibull(shape,scale) result(x)
        real(dp),intent(in)::shape,scale;real(dp)::u;call random_number(u);x=qinvweibull(u,shape,scale)
    end function rinvweibull
    pure real(dp) function minvweibull(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        if(order>=shape) then;m=huge(1.0_dp);else;m=scale**order*gamma(1.0_dp-order/shape);end if
    end function minvweibull

    pure real(dp) function dgumbel(x,alpha,scale) result(f)
        real(dp),intent(in)::x,alpha,scale;real(dp)::z
        if(scale<=0.0_dp) then;f=0.0_dp;else;z=(x-alpha)/scale;f=exp(-z-exp(-z))/scale;end if
    end function dgumbel
    pure real(dp) function pgumbel(x,alpha,scale) result(p)
        real(dp),intent(in)::x,alpha,scale;p=exp(-exp(-(x-alpha)/scale))
    end function pgumbel
    pure real(dp) function qgumbel(p,alpha,scale) result(x)
        real(dp),intent(in)::p,alpha,scale
        if(p<=0.0_dp) then;x=-huge(1.0_dp);else if(p>=1.0_dp) then;x=huge(1.0_dp);else;x=alpha-scale*log(-log(p));end if
    end function qgumbel
    real(dp) function rgumbel(alpha,scale) result(x)
        real(dp),intent(in)::alpha,scale;real(dp)::u;call random_number(u);x=qgumbel(u,alpha,scale)
    end function rgumbel
    pure real(dp) function mgumbel(order,alpha,scale) result(m)
        real(dp),intent(in)::order,alpha,scale;real(dp),parameter::euler=0.5772156649015328606_dp
        if(abs(order-1.0_dp)<1.0e-12_dp) then;m=alpha+euler*scale
        else if(abs(order-2.0_dp)<1.0e-12_dp) then;m=(pi*scale)**2/6.0_dp+(alpha+euler*scale)**2
        else;m=ieee_value(0.0_dp, ieee_quiet_nan);end if
    end function mgumbel
    pure real(dp) function mgfgumbel(t,alpha,scale) result(m)
        real(dp),intent(in)::t,alpha,scale
        if(scale*t>=1.0_dp) then;m=huge(1.0_dp);else;m=exp(alpha*t+log_gamma(1.0_dp-scale*t));end if
    end function mgfgumbel

    pure real(dp) function dlgamma(x,shapelog,ratelog) result(f)
        real(dp),intent(in)::x,shapelog,ratelog;real(dp)::y
        if(x<1.0_dp .or. min(shapelog,ratelog)<=0.0_dp) then;f=0.0_dp;return;end if
        y=log(x);f=ratelog**shapelog/gamma(shapelog)*y**(shapelog-1.0_dp)*exp(-ratelog*y)/x
    end function dlgamma
    pure real(dp) function plgamma(x,shapelog,ratelog) result(p)
        real(dp),intent(in)::x,shapelog,ratelog
        if(x<=1.0_dp) then;p=0.0_dp;else;p=reg_gamma_p(shapelog,ratelog*log(x));end if
    end function plgamma
    pure real(dp) function qlgamma(p,shapelog,ratelog) result(x)
        real(dp),intent(in)::p,shapelog,ratelog;integer::i;real(dp)::lo,hi,mid
        lo=1.0_dp;hi=exp(max(1.0_dp,shapelog/ratelog+10.0_dp*sqrt(shapelog)/ratelog))
        do i=1,100;mid=sqrt(lo*hi);if(plgamma(mid,shapelog,ratelog)<p) then;lo=mid;else;hi=mid;end if;end do;x=sqrt(lo*hi)
    end function qlgamma
    real(dp) function rlgamma(shapelog,ratelog) result(x)
        real(dp),intent(in)::shapelog,ratelog;x=exp(random_gamma(shapelog,1.0_dp/ratelog))
    end function rlgamma
    pure real(dp) function mlgamma(order,shapelog,ratelog) result(m)
        real(dp),intent(in)::order,shapelog,ratelog
        if(order>=ratelog) then;m=huge(1.0_dp);else;m=(1.0_dp-order/ratelog)**(-shapelog);end if
    end function mlgamma

    pure real(dp) function dtrgamma(x,shape1,shape2,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,scale;real(dp)::y
        if(x<=0.0_dp .or. min(shape1,shape2,scale)<=0.0_dp) then;f=0.0_dp;return;end if
        y=(x/scale)**shape2
        f=shape2/(x*gamma(shape1))*y**shape1*exp(-y)
    end function dtrgamma
    pure real(dp) function ptrgamma(x,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,scale
        if(x<=0.0_dp) then;p=0.0_dp;else;p=reg_gamma_p(shape1,(x/scale)**shape2);end if
    end function ptrgamma
    pure real(dp) function qtrgamma(p,shape1,shape2,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,scale;integer::i;real(dp)::lo,hi,mid
        lo=0.0_dp;hi=scale*max(2.0_dp,shape1**(1.0_dp/shape2)+10.0_dp)
        do while(ptrgamma(hi,shape1,shape2,scale)<p);hi=2.0_dp*hi;end do
        do i=1,100;mid=0.5_dp*(lo+hi);if(ptrgamma(mid,shape1,shape2,scale)<p) then;lo=mid;else;hi=mid;end if;end do;x=0.5_dp*(lo+hi)
    end function qtrgamma
    real(dp) function rtrgamma(shape1,shape2,scale) result(x)
        real(dp),intent(in)::shape1,shape2,scale;x=scale*random_gamma(shape1,1.0_dp)**(1.0_dp/shape2)
    end function rtrgamma
    pure real(dp) function mtrgamma(order,shape1,shape2,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,scale
        if(shape1+order/shape2<=0.0_dp) then;m=huge(1.0_dp);else;m=scale**order*gamma(shape1+order/shape2)/gamma(shape1);end if
    end function mtrgamma

    pure real(dp) function dtrbeta(x,shape1,shape2,shape3,scale) result(f)
        real(dp),intent(in)::x,shape1,shape2,shape3,scale
        real(dp)::v,u
        if(x<=0.0_dp .or. min(shape1,shape2,shape3,scale)<=0.0_dp) then
            f=0.0_dp
            return
        end if
        v=(x/scale)**shape2
        u=v/(1.0_dp+v)
        f=shape2*u**shape3*(1.0_dp-u)**shape1/(x*beta_fn(shape3,shape1))
    end function dtrbeta

    pure real(dp) function ptrbeta(x,shape1,shape2,shape3,scale) result(p)
        real(dp),intent(in)::x,shape1,shape2,shape3,scale
        real(dp)::v,u
        if(x<=0.0_dp) then
            p=0.0_dp
        else
            v=(x/scale)**shape2
            u=v/(1.0_dp+v)
            p=reg_beta(u,shape3,shape1)
        end if
    end function ptrbeta

    pure real(dp) function qtrbeta(p,shape1,shape2,shape3,scale) result(x)
        real(dp),intent(in)::p,shape1,shape2,shape3,scale
        real(dp)::u
        if(p<=0.0_dp) then
            x=0.0_dp
        else if(p>=1.0_dp) then
            x=huge(1.0_dp)
        else
            u=inv_reg_beta(p,shape3,shape1)
            x=scale*(u/(1.0_dp-u))**(1.0_dp/shape2)
        end if
    end function qtrbeta

    real(dp) function rtrbeta(shape1,shape2,shape3,scale) result(x)
        real(dp),intent(in)::shape1,shape2,shape3,scale
        real(dp)::u
        u=random_beta(shape3,shape1)
        x=scale*(u/(1.0_dp-u))**(1.0_dp/shape2)
    end function rtrbeta

    pure real(dp) function mtrbeta(order,shape1,shape2,shape3,scale) result(m)
        real(dp),intent(in)::order,shape1,shape2,shape3,scale
        real(dp)::t
        t=order/shape2
        if(t<=-shape3 .or. t>=shape1) then
            m=huge(1.0_dp)
        else
            m=scale**order*beta_fn(shape3+t,shape1-t)/beta_fn(shape1,shape3)
        end if
    end function mtrbeta

end module actuar_continuous
