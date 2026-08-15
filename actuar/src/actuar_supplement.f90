module actuar_supplement
    use actuar_kinds, only: dp
    use actuar_special, only: reg_gamma_p, reg_gamma_q, normal_cdf
    implicit none
    private
    public :: mexp_act, levexp_act, mgfexp_act
    public :: mgamma_act, levgamma_act, mgfgamma_act
    public :: mweibull_act, levweibull_act
    public :: mlnorm_act, levlnorm_act
    public :: mnorm_act, mgfnorm_act
    public :: munif_act, levunif_act, mgfunif_act
    public :: dinvgauss, pinvgauss, qinvgauss, rinvgauss, minvgauss, mgfinvgauss

contains

    pure real(dp) function mexp_act(order,rate) result(m)
        real(dp),intent(in)::order,rate
        if(rate<=0.0_dp .or. order<=-1.0_dp) then;m=huge(1.0_dp);else;m=gamma(order+1.0_dp)/rate**order;end if
    end function mexp_act

    pure real(dp) function levexp_act(limit,rate,order) result(m)
        real(dp),intent(in)::limit,rate,order
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        m=gamma(order+1.0_dp)/rate**order*reg_gamma_p(order+1.0_dp,rate*limit) &
          +limit**order*exp(-rate*limit)
    end function levexp_act

    pure real(dp) function mgfexp_act(t,rate) result(m)
        real(dp),intent(in)::t,rate
        if(t>=rate) then;m=huge(1.0_dp);else;m=rate/(rate-t);end if
    end function mgfexp_act

    pure real(dp) function mgamma_act(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        if(shape+order<=0.0_dp) then;m=huge(1.0_dp);else;m=scale**order*gamma(shape+order)/gamma(shape);end if
    end function mgamma_act

    pure real(dp) function levgamma_act(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        m=mgamma_act(order,shape,scale)*reg_gamma_p(shape+order,limit/scale) &
          +limit**order*reg_gamma_q(shape,limit/scale)
    end function levgamma_act

    pure real(dp) function mgfgamma_act(t,shape,scale) result(m)
        real(dp),intent(in)::t,shape,scale
        if(scale*t>=1.0_dp) then;m=huge(1.0_dp);else;m=(1.0_dp-scale*t)**(-shape);end if
    end function mgfgamma_act

    pure real(dp) function mweibull_act(order,shape,scale) result(m)
        real(dp),intent(in)::order,shape,scale
        if(1.0_dp+order/shape<=0.0_dp) then;m=huge(1.0_dp);else;m=scale**order*gamma(1.0_dp+order/shape);end if
    end function mweibull_act

    pure real(dp) function levweibull_act(limit,shape,scale,order) result(m)
        real(dp),intent(in)::limit,shape,scale,order
        real(dp)::z,a
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        z=(limit/scale)**shape;a=1.0_dp+order/shape
        m=scale**order*gamma(a)*reg_gamma_p(a,z)+limit**order*exp(-z)
    end function levweibull_act

    pure real(dp) function mlnorm_act(order,meanlog,sdlog) result(m)
        real(dp),intent(in)::order,meanlog,sdlog
        m=exp(order*meanlog+0.5_dp*(order*sdlog)**2)
    end function mlnorm_act

    pure real(dp) function levlnorm_act(limit,meanlog,sdlog,order) result(m)
        real(dp),intent(in)::limit,meanlog,sdlog,order
        real(dp)::z0,zr
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        z0=(log(limit)-meanlog)/sdlog
        zr=(log(limit)-meanlog-order*sdlog**2)/sdlog
        m=mlnorm_act(order,meanlog,sdlog)*normal_cdf(zr)+limit**order*(1.0_dp-normal_cdf(z0))
    end function levlnorm_act

    pure real(dp) function mnorm_act(order,mean,sd) result(m)
        integer,intent(in)::order
        real(dp),intent(in)::mean,sd
        real(dp)::m0,m1,m2
        integer::k
        if(order<0) then;m=huge(1.0_dp);return;end if
        if(order==0) then;m=1.0_dp;return;end if
        m0=1.0_dp;m1=mean
        if(order==1) then;m=m1;return;end if
        do k=2,order
            m2=mean*m1+real(k-1,dp)*sd**2*m0
            m0=m1;m1=m2
        end do
        m=m1
    end function mnorm_act

    pure real(dp) function mgfnorm_act(t,mean,sd) result(m)
        real(dp),intent(in)::t,mean,sd
        m=exp(mean*t+0.5_dp*sd**2*t**2)
    end function mgfnorm_act

    pure real(dp) function munif_act(order,xmin,xmax) result(m)
        real(dp),intent(in)::order,xmin,xmax
        if(abs(order+1.0_dp)<1.0e-14_dp) then
            m=log(xmax/xmin)/(xmax-xmin)
        else
            m=(xmax**(order+1.0_dp)-xmin**(order+1.0_dp))/((order+1.0_dp)*(xmax-xmin))
        end if
    end function munif_act

    pure real(dp) function levunif_act(limit,xmin,xmax,order) result(m)
        real(dp),intent(in)::limit,xmin,xmax,order
        real(dp)::u
        if(limit<=xmin) then;m=limit**order;return;end if
        if(limit>=xmax) then;m=munif_act(order,xmin,xmax);return;end if
        u=limit
        m=(u**(order+1.0_dp)-xmin**(order+1.0_dp))/((order+1.0_dp)*(xmax-xmin)) &
          +limit**order*(xmax-u)/(xmax-xmin)
    end function levunif_act

    pure real(dp) function mgfunif_act(t,xmin,xmax) result(m)
        real(dp),intent(in)::t,xmin,xmax
        if(abs(t)<1.0e-14_dp) then;m=1.0_dp;else;m=(exp(t*xmax)-exp(t*xmin))/(t*(xmax-xmin));end if
    end function mgfunif_act

    pure real(dp) function dinvgauss(x,mu,phi) result(f)
        real(dp),intent(in)::x,mu,phi
        real(dp)::u
        if(x<=0.0_dp .or. mu<=0.0_dp .or. phi<=0.0_dp) then;f=0.0_dp;return;end if
        u=(x-mu)/mu
        f=exp(-u*u/(2.0_dp*phi*x))/sqrt(2.0_dp*acos(-1.0_dp)*phi*x**3)
    end function dinvgauss

    pure real(dp) function pinvgauss(x,mu,phi) result(p)
        real(dp),intent(in)::x,mu,phi
        real(dp)::r,a,b,phim,qm
        if(x<=0.0_dp) then;p=0.0_dp;return;end if
        qm=x/mu;phim=phi*mu;r=sqrt(x*phi)
        a=(qm-1.0_dp)/r;b=-(qm+1.0_dp)/r
        p=normal_cdf(a)+exp(2.0_dp/phim)*normal_cdf(b)
        p=max(0.0_dp,min(1.0_dp,p))
    end function pinvgauss

    pure real(dp) function qinvgauss(p,mu,phi) result(x)
        real(dp),intent(in)::p,mu,phi
        real(dp)::lo,hi,mid
        integer::i
        if(p<=0.0_dp) then;x=0.0_dp;return;end if
        if(p>=1.0_dp) then;x=huge(1.0_dp);return;end if
        lo=0.0_dp;hi=max(mu,1.0_dp)
        do while(pinvgauss(hi,mu,phi)<p);hi=2.0_dp*hi;end do
        do i=1,110
            mid=0.5_dp*(lo+hi)
            if(pinvgauss(mid,mu,phi)<p) then;lo=mid;else;hi=mid;end if
        end do
        x=0.5_dp*(lo+hi)
    end function qinvgauss

    real(dp) function rinvgauss(mu,phi) result(x)
        real(dp),intent(in)::mu,phi
        real(dp)::v,y,u,trial,z,u1,u2
        call random_number(u1);call random_number(u2)
        z=sqrt(-2.0_dp*log(max(u1,tiny(1.0_dp))))*cos(2.0_dp*acos(-1.0_dp)*u2)
        y=z*z;v=mu*phi
        trial=mu+0.5_dp*mu*v*y-0.5_dp*mu*sqrt(4.0_dp*v*y+v*v*y*y)
        call random_number(u)
        if(u<=mu/(mu+trial)) then;x=trial;else;x=mu*mu/trial;end if
    end function rinvgauss

    pure real(dp) function minvgauss(order,mu,phi) result(m)
        integer,intent(in)::order
        real(dp),intent(in)::mu,phi
        integer::k
        real(dp)::term,sumv
        if(order<0) then;m=huge(1.0_dp);return;end if
        if(order==0) then;m=1.0_dp;return;end if
        sumv=0.0_dp
        do k=0,order-1
            term=exp(log_gamma(real(order+k,dp))-log_gamma(real(k+1,dp)) &
                 -log_gamma(real(order-k,dp)))*(0.5_dp*mu*phi)**k
            sumv=sumv+term
        end do
        m=mu**order*sumv
    end function minvgauss

    pure real(dp) function mgfinvgauss(t,mu,phi) result(m)
        real(dp),intent(in)::t,mu,phi
        real(dp)::arg
        arg=1.0_dp-2.0_dp*phi*mu*mu*t
        if(arg<0.0_dp) then;m=huge(1.0_dp);else;m=exp((1.0_dp-sqrt(arg))/(phi*mu));end if
    end function mgfinvgauss

end module actuar_supplement
