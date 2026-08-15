module actuar_supplement_v02
    use actuar_kinds, only: dp
    use actuar_special, only: beta_fn, reg_beta, reg_gamma_p, reg_gamma_q, normal_cdf
    use actuar_supplement, only: pinvgauss
    use expint_mod, only: gamma_inc
    implicit none
    private
    public :: mbeta_act, levbeta_act
    public :: mchisq_act, levchisq_act, mgfchisq_act
    public :: levinvexp_act, levinvgauss_act
contains
    pure real(dp) function mbeta_act(order,shape1,shape2) result(m)
        real(dp),intent(in)::order,shape1,shape2
        if(order<=-shape1) then
            m=huge(1.0_dp)
        else
            m=beta_fn(shape1+order,shape2)/beta_fn(shape1,shape2)
        end if
    end function mbeta_act

    pure real(dp) function levbeta_act(limit,shape1,shape2,order) result(m)
        real(dp),intent(in)::limit,shape1,shape2,order
        real(dp)::tmp,l
        if(order<=-shape1) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        if(limit>=1.0_dp) then;m=mbeta_act(order,shape1,shape2);return;end if
        tmp=shape1+order;l=limit
        m=beta_fn(tmp,shape2)/beta_fn(shape1,shape2)*reg_beta(l,tmp,shape2) &
          + limit**order*(1.0_dp-reg_beta(l,shape1,shape2))
    end function levbeta_act

    pure real(dp) function mchisq_act(order,df,ncp) result(m)
        real(dp),intent(in)::order,df,ncp
        real(dp),allocatable::res(:)
        integer::i,j,n
        if(order<=-df/2.0_dp) then;m=huge(1.0_dp);return;end if
        if(order==0.0_dp) then;m=1.0_dp;return;end if
        if(ncp==0.0_dp) then
            m=2.0_dp**order*gamma(order+df/2.0_dp)/gamma(df/2.0_dp);return
        end if
        if(order<1.0_dp .or. abs(order-real(nint(order),dp))>64.0_dp*epsilon(1.0_dp)) then
            m=huge(1.0_dp);return
        end if
        n=nint(order);allocate(res(0:n));res=0.0_dp;res(0)=1.0_dp;res(1)=df+ncp
        do i=2,n
            res(i)=2.0_dp**(i-1)*(df+real(i,dp)*ncp)
            do j=1,i-1
                res(i)=res(i)+2.0_dp**(j-1)*(df+real(j,dp)*ncp)*res(i-j)/gamma(real(i-j+1,dp))
            end do
            res(i)=res(i)*gamma(real(i,dp))
        end do
        m=res(n)
    end function mchisq_act

    pure real(dp) function levchisq_act(limit,df,ncp,order) result(m)
        real(dp),intent(in)::limit,df,ncp,order
        real(dp)::tmp,u
        if(order<=-df/2.0_dp) then;m=huge(1.0_dp);return;end if
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        if(ncp/=0.0_dp) then;m=huge(1.0_dp);return;end if
        tmp=order+df/2.0_dp;u=limit/2.0_dp
        m=2.0_dp**order*gamma(tmp)*reg_gamma_p(tmp,u)/gamma(df/2.0_dp) &
          + limit**order*reg_gamma_q(df/2.0_dp,u)
    end function levchisq_act

    pure real(dp) function mgfchisq_act(t,df,ncp) result(mgf)
        real(dp),intent(in)::t,df,ncp
        if(2.0_dp*t>=1.0_dp) then;mgf=huge(1.0_dp);return;end if
        mgf=exp(ncp*t/(1.0_dp-2.0_dp*t)-0.5_dp*df*log(1.0_dp-2.0_dp*t))
    end function mgfchisq_act

    pure real(dp) function levinvexp_act(limit,scale,order) result(m)
        real(dp),intent(in)::limit,scale,order
        real(dp)::u
        if(limit<=0.0_dp) then;m=0.0_dp;return;end if
        u=scale/limit
        m=scale**order*gamma_inc(1.0_dp-order,u)+limit**order*(1.0_dp-exp(-u))
    end function levinvexp_act

    pure real(dp) function levinvgauss_act(limit,mu,phi) result(m)
        real(dp),intent(in)::limit,mu,phi
        real(dp)::xm,phim,r,x,z2,a,b,tail
        if(limit<=0.0_dp .or. phi<=0.0_dp) then;m=0.0_dp;return;end if
        xm=limit/mu;phim=phi*mu;r=sqrt(limit*phi)
        x=(xm-1.0_dp)/r;z2=-(xm+1.0_dp)/r
        a=normal_cdf(x);b=exp(2.0_dp/phim)*normal_cdf(z2)
        tail=1.0_dp-pinvgauss(limit,mu,phi)
        m=mu*max(0.0_dp,a-b)+limit*max(0.0_dp,tail)
    end function levinvgauss_act
end module actuar_supplement_v02
