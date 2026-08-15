module actuar_continuous_v02
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use actuar_kinds, only: dp
    use actuar_special, only: beta_fn, reg_beta, inv_reg_beta, reg_gamma_p, reg_gamma_q, random_beta, random_gamma
    use actuar_special_v02, only: betaint_raw, choose_int
    use actuar_continuous, only: dpareto, ppareto, mpareto, levpareto, dllogis, pllogis, mllogis, &
                                  dburr, pburr, mburr, dinvburr, pinvburr, minvburr, &
                                  dgenpareto, pgenpareto, dgenbeta, pgenbeta, &
                                  dinvgamma, pinvgamma, dinvweibull, pinvweibull, &
                                  dtrgamma, ptrgamma, mtrgamma, dtrbeta, ptrbeta, mtrbeta
    use expint_mod, only: gamma_inc
    implicit none
    private
    public :: dfpareto, pfpareto, qfpareto, rfpareto, mfpareto, levfpareto
    public :: dinvpareto, pinvpareto, qinvpareto, rinvpareto, minvpareto, levinvpareto
    public :: dinvtrgamma, pinvtrgamma, qinvtrgamma, rinvtrgamma, minvtrgamma, levinvtrgamma
    public :: mpareto2, levpareto2, mpareto3, levpareto3, mpareto4, levpareto4
    public :: levburr, levinvburr, levllogis
    public :: mparalogis, levparalogis, minvparalogis, levinvparalogis
    public :: levgenpareto, levgenbeta, levtrgamma, levtrbeta
    public :: levinvgamma, levinvweibull
contains
    pure real(dp) function quiet_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function quiet_nan

    pure real(dp) function safe_pow_limit(limit, order) result(v)
        real(dp), intent(in) :: limit, order
        if (limit <= 0.0_dp .and. order == 0.0_dp) then
            v = 1.0_dp
        else
            v = limit**order
        end if
    end function safe_pow_limit

    pure real(dp) function dfpareto(x, xmin, shape1, shape2, shape3, scale) result(f)
        real(dp), intent(in) :: x, xmin, shape1, shape2, shape3, scale
        real(dp) :: v, u
        if (shape1 <= 0.0_dp .or. shape2 <= 0.0_dp .or. shape3 <= 0.0_dp .or. scale <= 0.0_dp) then
            f = quiet_nan(); return
        end if
        if (x < xmin) then
            f = 0.0_dp; return
        end if
        if (x == xmin) then
            if (shape2*shape3 < 1.0_dp) then
                f = huge(1.0_dp)
            else if (shape2*shape3 > 1.0_dp) then
                f = 0.0_dp
            else
                f = shape2/(scale*beta_fn(shape3,shape1))
            end if
            return
        end if
        v = ((x-xmin)/scale)**shape2
        u = v/(1.0_dp+v)
        f = shape2*u**shape3*(1.0_dp-u)**shape1/((x-xmin)*beta_fn(shape3,shape1))
    end function dfpareto

    pure real(dp) function pfpareto(x, xmin, shape1, shape2, shape3, scale) result(p)
        real(dp), intent(in) :: x, xmin, shape1, shape2, shape3, scale
        real(dp) :: v, u
        if (x <= xmin) then
            p = 0.0_dp
            return
        end if
        v = ((x-xmin)/scale)**shape2
        u = v/(1.0_dp+v)
        p = reg_beta(u,shape3,shape1)
    end function pfpareto

    pure real(dp) function qfpareto(p, xmin, shape1, shape2, shape3, scale) result(x)
        real(dp), intent(in) :: p, xmin, shape1, shape2, shape3, scale
        real(dp) :: u
        if (p <= 0.0_dp) then
            x = xmin
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp)
        else
            u = inv_reg_beta(p,shape3,shape1)
            x = xmin + scale*(u/(1.0_dp-u))**(1.0_dp/shape2)
        end if
    end function qfpareto

    real(dp) function rfpareto(xmin, shape1, shape2, shape3, scale) result(x)
        real(dp), intent(in) :: xmin, shape1, shape2, shape3, scale
        real(dp) :: u
        u = random_beta(shape3,shape1)
        x = xmin + scale*(u/(1.0_dp-u))**(1.0_dp/shape2)
    end function rfpareto

    pure real(dp) function mfpareto(order, xmin, shape1, shape2, shape3, scale) result(m)
        real(dp), intent(in) :: order, xmin, shape1, shape2, shape3, scale
        integer :: i, n
        real(dp) :: tmp, sumv, r, be
        if (xmin == 0.0_dp) then
            m = mtrbeta(order,shape1,shape2,shape3,scale)
            return
        end if
        if (order < 0.0_dp) then
            m = quiet_nan(); return
        end if
        if (order >= shape1*shape2) then
            m = huge(1.0_dp); return
        end if
        n = nint(order)
        be = beta_fn(shape1,shape3)
        r = scale/xmin
        sumv = be
        do i = 1, n
            tmp = real(i,dp)/shape2
            sumv = sumv + choose_int(n,i)*r**i*beta_fn(shape3+tmp,shape1-tmp)
        end do
        m = xmin**n*sumv/be
    end function mfpareto

    real(dp) function levfpareto(limit, xmin, shape1, shape2, shape3, scale, order) result(m)
        real(dp), intent(in) :: limit, xmin, shape1, shape2, shape3, scale, order
        integer :: i, n
        real(dp) :: u, tmp, sumv, r
        if (limit <= xmin) then
            m = 0.0_dp; return
        end if
        if (xmin == 0.0_dp) then
            m = levtrbeta(limit,shape1,shape2,shape3,scale,order)
            return
        end if
        if (order < 0.0_dp) then
            m = quiet_nan(); return
        end if
        n = nint(order)
        u = ((limit-xmin)/scale)**shape2
        u = u/(1.0_dp+u)
        r = scale/xmin
        sumv = betaint_raw(u,shape3,shape1)
        do i = 1, n
            tmp = real(i,dp)/shape2
            sumv = sumv + choose_int(n,i)*r**i*betaint_raw(u,shape3+tmp,shape1-tmp)
        end do
        m = xmin**n*sumv/(gamma(shape1)*gamma(shape3)) + limit**n*(1.0_dp-pfpareto(limit,xmin,shape1,shape2,shape3,scale))
    end function levfpareto

    pure real(dp) function dinvpareto(x, shape, scale) result(f)
        real(dp), intent(in) :: x, shape, scale
        real(dp) :: u
        if (shape <= 0.0_dp .or. scale <= 0.0_dp) then
            f = quiet_nan(); return
        end if
        if (x < 0.0_dp) then
            f = 0.0_dp; return
        end if
        if (x == 0.0_dp) then
            if (shape < 1.0_dp) then
                f = huge(1.0_dp)
            else if (shape > 1.0_dp) then
                f = 0.0_dp
            else
                f = 1.0_dp/scale
            end if
            return
        end if
        u = x/(x+scale)
        f = shape*u**shape*(1.0_dp-u)/x
    end function dinvpareto

    pure real(dp) function pinvpareto(x, shape, scale) result(p)
        real(dp), intent(in) :: x, shape, scale
        if (x <= 0.0_dp) then
            p = 0.0_dp
        else
            p = (x/(x+scale))**shape
        end if
    end function pinvpareto

    pure real(dp) function qinvpareto(p, shape, scale) result(x)
        real(dp), intent(in) :: p, shape, scale
        if (p <= 0.0_dp) then
            x = 0.0_dp
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp)
        else
            x = scale/(p**(-1.0_dp/shape)-1.0_dp)
        end if
    end function qinvpareto

    real(dp) function rinvpareto(shape, scale) result(x)
        real(dp), intent(in) :: shape, scale
        real(dp) :: u
        call random_number(u)
        x = qinvpareto(u,shape,scale)
    end function rinvpareto

    pure real(dp) function minvpareto(order, shape, scale) result(m)
        real(dp), intent(in) :: order, shape, scale
        if (order <= -shape .or. order >= 1.0_dp) then
            m = huge(1.0_dp)
        else
            m = scale**order*gamma(shape+order)*gamma(1.0_dp-order)/gamma(shape)
        end if
    end function minvpareto

    real(dp) function levinvpareto(limit, shape, scale, order) result(m)
        real(dp), intent(in) :: limit, shape, scale, order
        real(dp) :: u, a, b, integ
        if (order <= -shape) then
            m = huge(1.0_dp); return
        end if
        if (limit <= 0.0_dp) then
            m = 0.0_dp; return
        end if
        u = limit/(limit+scale)
        a = shape+order
        b = 1.0_dp-order
        if (b > 0.0_dp .or. abs(b-real(nint(b),dp)) > 64.0_dp*epsilon(1.0_dp)) then
            integ = betaint_raw(u,a,b)/gamma(shape+1.0_dp)
        else
            integ = invpareto_integral(u,shape,order)
        end if
        m = scale**order*shape*integ + limit**order*(1.0_dp-pinvpareto(limit,shape,scale))
    end function levinvpareto

    real(dp) function invpareto_integral(u, shape, order) result(v)
        real(dp), intent(in) :: u, shape, order
        integer, parameter :: n = 20000
        integer :: i
        real(dp) :: h, x, s
        h = u/real(n,dp)
        s = 0.0_dp
        do i = 1, n
            x = (real(i,dp)-0.5_dp)*h
            s = s + x**(shape+order-1.0_dp)*(1.0_dp-x)**(-order)
        end do
        v = h*s
    end function invpareto_integral

    pure real(dp) function dinvtrgamma(x, shape1, shape2, scale) result(f)
        real(dp), intent(in) :: x, shape1, shape2, scale
        real(dp) :: u
        if (x <= 0.0_dp) then
            f = 0.0_dp; return
        end if
        u = (scale/x)**shape2
        f = shape2*u**shape1*exp(-u)/(x*gamma(shape1))
    end function dinvtrgamma

    pure real(dp) function pinvtrgamma(x, shape1, shape2, scale) result(p)
        real(dp), intent(in) :: x, shape1, shape2, scale
        real(dp) :: u
        if (x <= 0.0_dp) then
            p = 0.0_dp
        else
            u = (scale/x)**shape2
            p = reg_gamma_q(shape1,u)
        end if
    end function pinvtrgamma

    pure real(dp) function qinvtrgamma(p, shape1, shape2, scale) result(x)
        real(dp), intent(in) :: p, shape1, shape2, scale
        real(dp) :: lo, hi, mid
        integer :: i
        if (p <= 0.0_dp) then
            x = 0.0_dp; return
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp); return
        end if
        lo = 0.0_dp
        hi = max(scale,1.0_dp)
        do while (pinvtrgamma(hi,shape1,shape2,scale) < p)
            hi = 2.0_dp*hi
            if (hi > huge(1.0_dp)/4.0_dp) exit
        end do
        do i = 1, 120
            mid = 0.5_dp*(lo+hi)
            if (pinvtrgamma(mid,shape1,shape2,scale) < p) then
                lo = mid
            else
                hi = mid
            end if
        end do
        x = 0.5_dp*(lo+hi)
    end function qinvtrgamma

    real(dp) function rinvtrgamma(shape1, shape2, scale) result(x)
        real(dp), intent(in) :: shape1, shape2, scale
        x = scale*random_gamma(shape1,1.0_dp)**(-1.0_dp/shape2)
    end function rinvtrgamma

    pure real(dp) function minvtrgamma(order, shape1, shape2, scale) result(m)
        real(dp), intent(in) :: order, shape1, shape2, scale
        if (order >= shape1*shape2) then
            m = huge(1.0_dp)
        else
            m = scale**order*gamma(shape1-order/shape2)/gamma(shape1)
        end if
    end function minvtrgamma

    pure real(dp) function levinvtrgamma(limit, shape1, shape2, scale, order) result(m)
        real(dp), intent(in) :: limit, shape1, shape2, scale, order
        real(dp) :: u
        if (limit <= 0.0_dp) then
            m = 0.0_dp; return
        end if
        u = (scale/limit)**shape2
        m = scale**order*gamma_inc(shape1-order/shape2,u)/gamma(shape1) &
            + limit**order*reg_gamma_p(shape1,u)
    end function levinvtrgamma

    pure real(dp) function mpareto2(order, xmin, shape, scale) result(m)
        real(dp), intent(in) :: order, xmin, shape, scale
        integer :: i, n
        real(dp) :: sumv, r
        if (xmin == 0.0_dp) then
            m = mpareto(order,shape,scale); return
        end if
        if (order < 0.0_dp) then
            m = quiet_nan(); return
        end if
        if (order >= shape) then
            m = huge(1.0_dp); return
        end if
        n = nint(order)
        r = scale/xmin
        sumv = gamma(shape)
        do i=1,n
            sumv = sumv + choose_int(n,i)*r**i*gamma(1.0_dp+real(i,dp))*gamma(shape-real(i,dp))
        end do
        m = xmin**n*sumv/gamma(shape)
    end function mpareto2

    real(dp) function levpareto2(limit, xmin, shape, scale, order) result(m)
        real(dp), intent(in) :: limit, xmin, shape, scale, order
        integer :: i, n
        real(dp) :: u, sumv, r
        if (limit <= xmin) then
            m = 0.0_dp; return
        end if
        if (xmin == 0.0_dp) then
            m = levpareto(limit,shape,scale,order); return
        end if
        if (order < 0.0_dp) then
            m = quiet_nan(); return
        end if
        n = nint(order)
        u = scale/(limit-xmin+scale)
        r = scale/xmin
        sumv = betaint_raw(1.0_dp-u,1.0_dp,shape)
        do i=1,n
            sumv = sumv + choose_int(n,i)*r**i*betaint_raw(1.0_dp-u,1.0_dp+real(i,dp),shape-real(i,dp))
        end do
        m = xmin**n*sumv/gamma(shape) + limit**n*(1.0_dp-pareto2_cdf(limit,xmin,shape,scale))
    end function levpareto2

    pure real(dp) function pareto2_cdf(x,xmin,shape,scale) result(p)
        real(dp), intent(in) :: x,xmin,shape,scale
        if (x <= xmin) then
            p=0.0_dp
        else
            p=ppareto(x-xmin,shape,scale)
        end if
    end function pareto2_cdf

    pure real(dp) function mpareto3(order, xmin, shape, scale) result(m)
        real(dp), intent(in) :: order, xmin, shape, scale
        integer :: i, n
        real(dp) :: tmp, sumv, r
        if (xmin == 0.0_dp) then
            m = mllogis(order,shape,scale); return
        end if
        if (order < 0.0_dp) then
            m = quiet_nan(); return
        end if
        if (order >= shape) then
            m = huge(1.0_dp); return
        end if
        n=nint(order); r=scale/xmin; sumv=1.0_dp
        do i=1,n
            tmp=real(i,dp)/shape
            sumv=sumv+choose_int(n,i)*r**i*gamma(1.0_dp+tmp)*gamma(1.0_dp-tmp)
        end do
        m=xmin**n*sumv
    end function mpareto3

    real(dp) function levpareto3(limit, xmin, shape, scale, order) result(m)
        real(dp), intent(in) :: limit, xmin, shape, scale, order
        integer :: i,n
        real(dp) :: v,u,tmp,sumv,r
        if(limit<=xmin) then;m=0.0_dp;return;end if
        if(xmin==0.0_dp) then;m=levllogis(limit,shape,scale,order);return;end if
        if(order<0.0_dp) then;m=quiet_nan();return;end if
        n=nint(order);r=scale/xmin
        v=((limit-xmin)/scale)**shape;u=v/(1.0_dp+v)
        sumv=betaint_raw(u,1.0_dp,1.0_dp)
        do i=1,n
            tmp=real(i,dp)/shape
            sumv=sumv+choose_int(n,i)*r**i*betaint_raw(u,1.0_dp+tmp,1.0_dp-tmp)
        end do
        m=xmin**n*sumv+limit**n*(1.0_dp-loglogis_shift_cdf(limit,xmin,shape,scale))
    end function levpareto3

    pure real(dp) function loglogis_shift_cdf(x,xmin,shape,scale) result(p)
        real(dp),intent(in)::x,xmin,shape,scale
        if(x<=xmin) then;p=0.0_dp;else;p=pllogis(x-xmin,shape,scale);end if
    end function loglogis_shift_cdf

    pure real(dp) function mpareto4(order, xmin, shape1, shape2, scale) result(m)
        real(dp), intent(in) :: order, xmin, shape1, shape2, scale
        integer :: i,n
        real(dp)::tmp,sumv,r
        if(xmin==0.0_dp) then;m=mburr(order,shape1,shape2,scale);return;end if
        if(order<0.0_dp) then;m=quiet_nan();return;end if
        if(order>=shape1*shape2) then;m=huge(1.0_dp);return;end if
        n=nint(order);r=scale/xmin;sumv=gamma(shape1)
        do i=1,n
            tmp=real(i,dp)/shape2
            sumv=sumv+choose_int(n,i)*r**i*gamma(1.0_dp+tmp)*gamma(shape1-tmp)
        end do
        m=xmin**n*sumv/gamma(shape1)
    end function mpareto4

    real(dp) function levpareto4(limit, xmin, shape1, shape2, scale, order) result(m)
        real(dp), intent(in) :: limit, xmin, shape1, shape2, scale, order
        integer::i,n
        real(dp)::v,u,tmp,sumv,r
        if(limit<=xmin) then;m=0.0_dp;return;end if
        if(xmin==0.0_dp) then;m=levburr(limit,shape1,shape2,scale,order);return;end if
        if(order<0.0_dp) then;m=quiet_nan();return;end if
        n=nint(order);r=scale/xmin
        v=((limit-xmin)/scale)**shape2;u=1.0_dp/(1.0_dp+v)
        sumv=betaint_raw(1.0_dp-u,1.0_dp,shape1)
        do i=1,n
            tmp=real(i,dp)/shape2
            sumv=sumv+choose_int(n,i)*r**i*betaint_raw(1.0_dp-u,1.0_dp+tmp,shape1-tmp)
        end do
        m=xmin**n*sumv/gamma(shape1)+limit**n*(1.0_dp-burr_shift_cdf(limit,xmin,shape1,shape2,scale))
    end function levpareto4

    pure real(dp) function burr_shift_cdf(x,xmin,shape1,shape2,scale) result(p)
        real(dp),intent(in)::x,xmin,shape1,shape2,scale
        if(x<=xmin) then;p=0.0_dp;else;p=pburr(x-xmin,shape1,shape2,scale);end if
    end function burr_shift_cdf

    pure real(dp) function levburr(limit,shape1,shape2,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,scale,order
        real(dp)::v,u,tmp
        if(order<=-shape2) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape2;v=(limit/scale)**shape2;u=v/(1.0_dp+v)
        m=scale**order*betaint_raw(u,1.0_dp+tmp,shape1-tmp)/gamma(shape1) &
          + limit**order*(1.0_dp-pburr(limit,shape1,shape2,scale))
    end function levburr

    pure real(dp) function levinvburr(limit,shape1,shape2,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,scale,order
        real(dp)::v,u,tmp
        if(order<=-shape1*shape2) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape2;v=(limit/scale)**shape2;u=v/(1.0_dp+v)
        m=scale**order*betaint_raw(u,shape1+tmp,1.0_dp-tmp)/gamma(shape1) &
          + limit**order*(1.0_dp-pinvburr(limit,shape1,shape2,scale))
    end function levinvburr

    pure real(dp) function levllogis(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::v,u,tmp
        if(order<=-shape) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape;v=(limit/scale)**shape;u=v/(1.0_dp+v)
        m=scale**order*betaint_raw(u,1.0_dp+tmp,1.0_dp-tmp) &
          + limit**order*(1.0_dp-pllogis(limit,shape,scale))
    end function levllogis

    pure real(dp) function mparalogis(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        real(dp)::tmp
        if(order<=-shape .or. order>=shape*shape) then;m=huge(1.0_dp);return;end if
        tmp=order/shape
        m=scale**order*gamma(1.0_dp+tmp)*gamma(shape-tmp)/gamma(shape)
    end function mparalogis

    pure real(dp) function levparalogis(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::v,u,tmp
        if(order<=-shape) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape;v=(limit/scale)**shape;u=v/(1.0_dp+v)
        m=scale**order*betaint_raw(u,1.0_dp+tmp,shape-tmp)/gamma(shape) &
          + limit**order*(1.0_dp-paralogis_cdf(limit,shape,scale))
    end function levparalogis

    pure real(dp) function paralogis_cdf(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale
        p=1.0_dp-(1.0_dp+(x/scale)**shape)**(-shape)
    end function paralogis_cdf

    pure real(dp) function minvparalogis(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        real(dp)::tmp
        if(order<=-shape*shape .or. order>=shape) then;m=huge(1.0_dp);return;end if
        tmp=order/shape
        m=scale**order*gamma(shape+tmp)*gamma(1.0_dp-tmp)/gamma(shape)
    end function minvparalogis

    pure real(dp) function levinvparalogis(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::v,u,tmp,p
        if(order<=-shape*shape) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape;v=(limit/scale)**shape;u=v/(1.0_dp+v)
        p=u**shape
        m=scale**order*betaint_raw(u,shape+tmp,1.0_dp-tmp)/gamma(shape) &
          + limit**order*(1.0_dp-p)
    end function levinvparalogis

    pure real(dp) function levgenpareto(limit,shape1,shape2,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,scale,order
        real(dp)::u
        if(order<=-shape2) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        u=limit/(limit+scale)
        m=scale**order*betaint_raw(u,shape2+order,shape1-order)/(gamma(shape1)*gamma(shape2)) &
          + limit**order*(1.0_dp-pgenpareto(limit,shape1,shape2,scale))
    end function levgenpareto

    pure real(dp) function levgenbeta(limit,shape1,shape2,shape3,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,shape3,scale,order
        real(dp)::u,tmp
        if(order<=-shape1*shape3) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        if(limit>=scale) then
            m=scale**order*beta_fn(shape1+order/shape3,shape2)/beta_fn(shape1,shape2)
            return
        end if
        tmp=order/shape3;u=(limit/scale)**shape3
        m=scale**order*beta_fn(shape1+tmp,shape2)/beta_fn(shape1,shape2)*reg_beta(u,shape1+tmp,shape2) &
          + limit**order*(1.0_dp-pgenbeta(limit,shape1,shape2,shape3,scale))
    end function levgenbeta

    pure real(dp) function levtrgamma(limit,shape1,shape2,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,scale,order
        real(dp)::u,tmp
        if(order<=-shape1*shape2) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=shape1+order/shape2;u=(limit/scale)**shape2
        m=scale**order*gamma(tmp)*reg_gamma_p(tmp,u)/gamma(shape1) &
          + limit**order*(1.0_dp-ptrgamma(limit,shape1,shape2,scale))
    end function levtrgamma

    pure real(dp) function levtrbeta(limit,shape1,shape2,shape3,scale,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,shape3,scale,order
        real(dp)::v,u,tmp
        if(order<=-shape3*shape2) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        tmp=order/shape2;v=(limit/scale)**shape2;u=v/(1.0_dp+v)
        m=scale**order*betaint_raw(u,shape3+tmp,shape1-tmp)/(gamma(shape1)*gamma(shape3)) &
          + limit**order*(1.0_dp-ptrbeta(limit,shape1,shape2,shape3,scale))
    end function levtrbeta

    pure real(dp) function levinvgamma(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::u
        if(order>=shape) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        u=scale/limit
        m=scale**order*gamma_inc(shape-order,u)/gamma(shape) + limit**order*(1.0_dp-pingamma_cdf(limit,shape,scale))
    end function levinvgamma

    pure real(dp) function pingamma_cdf(x,shape,scale) result(p)
        real(dp),intent(in)::x,shape,scale
        p=pinvgamma(x,shape,scale)
    end function pingamma_cdf

    pure real(dp) function levinvweibull(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::u
        if(order>=shape) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        u=(scale/limit)**shape
        m=scale**order*gamma_inc(1.0_dp-order/shape,u) + limit**order*(1.0_dp-pinvweibull(limit,shape,scale))
    end function levinvweibull
end module actuar_continuous_v02
