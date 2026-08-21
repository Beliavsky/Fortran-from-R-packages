module boot_profiles
    use boot_kinds, only : dp
    implicit none
    private
    public :: empirical_loglikelihood, eef_loglikelihood, likelihood_ci
contains
    subroutine empirical_loglikelihood(u,loglik,lambda,info)
        real(dp),intent(in)::u(:)
        real(dp),intent(out)::loglik,lambda
        integer,intent(out),optional::info
        real(dp)::lo,hi,mid,f,fl,fh,den
        integer::i,it
        lo=-huge(1.0_dp)
        hi=huge(1.0_dp)
        do i=1,size(u)
            if(u(i)>0.0_dp)lo=max(lo,-1.0_dp/u(i))
            if(u(i)<0.0_dp)hi=min(hi,-1.0_dp/u(i))
        end do
        if(lo<=-huge(1.0_dp)/2)lo=-1.0e6_dp
        if(hi>=huge(1.0_dp)/2)hi=1.0e6_dp
        lo=lo+1.0e-12_dp*max(1.0_dp,abs(lo))
        hi=hi-1.0e-12_dp*max(1.0_dp,abs(hi))
        fl=sum(u/(1.0_dp+lo*u))
        fh=sum(u/(1.0_dp+hi*u))
        if(fl*fh>0.0_dp)then
            lambda=0.0_dp
            loglik=-huge(1.0_dp)
            if(present(info))info=1
            return
        end if
        do it=1,200
            mid=0.5_dp*(lo+hi)
            f=sum(u/(1.0_dp+mid*u))
            if(abs(f)<1.0e-12_dp)exit
            if(fl*f<=0.0_dp)then
            hi=mid
            fh=f
            else
            lo=mid
            fl=f
            end if
        end do
        lambda=mid
        den=1.0_dp+lambda*minval(u)
        if(den<=0.0_dp)then
        loglik=-huge(1.0_dp)
        if(present(info))info=2
        return
        end if
        loglik=-sum(log(1.0_dp+lambda*u))
        if(present(info))info=0
    end subroutine empirical_loglikelihood

    subroutine eef_loglikelihood(u,loglik,eef,lambda,info)
        real(dp),intent(in)::u(:)
        real(dp),intent(out)::loglik,eef,lambda
        integer,intent(out),optional::info
        real(dp)::lam,f,der,mx,sf
        real(dp),allocatable::fit(:)
        integer::it
        allocate(fit(size(u)))
        lam=0.0_dp
        do it=1,100
            mx=maxval(lam*u)
            fit=exp(lam*u-mx)
            f=sum(u*fit)
            der=sum(u*u*fit)
            if(abs(f)<=1.0e-12_dp*max(1.0_dp,der))exit
            if(der<=tiny(1.0_dp))exit
            lam=lam-f/der
        end do
        lambda=lam
        fit=exp(lam*u)
        sf=sum(fit)
        loglik=sum(log(fit)-log(sf))
        eef=sum(fit-1.0_dp)
        if(present(info))info=merge(0,1,it<=100)
    end subroutine eef_loglikelihood

    subroutine likelihood_ci(theta,like,limit,lo,hi)
        real(dp),intent(in)::theta(:),like(:),limit
        real(dp),intent(out)::lo,hi
        integer::i,n,i0,j0
        real(dp)::x(3),y(3),a,b,c,disc
        n=size(theta)
        if(size(like)/=n)error stop "likelihood_ci: mismatch"
        i0=0
        j0=0
        do i=1,n
        if(like(i)>limit)then
        if(i0==0)i0=i
        j0=i
        end if
        end do
        if(i0==0 .or. i0==j0)error stop "likelihood_ci: insufficient crossing"
        if(i0==1)then
        lo=-huge(1.0_dp)
        else
            x=theta(i0-1:i0+1)
            y=like(i0-1:i0+1)-limit
            call quadratic_fit(x,y,a,b,c)
            disc=b*b-4*a*c
            lo=(-b+sqrt(max(0.0_dp,disc)))/(2*a)
        end if
        if(j0==n)then
        hi=huge(1.0_dp)
        else
            x=theta(j0-1:j0+1)
            y=like(j0-1:j0+1)-limit
            call quadratic_fit(x,y,a,b,c)
            disc=b*b-4*a*c
            hi=(-b-sqrt(max(0.0_dp,disc)))/(2*a)
        end if
    end subroutine likelihood_ci

    subroutine quadratic_fit(x,y,a,b,c)
        real(dp),intent(in)::x(3),y(3)
        real(dp),intent(out)::a,b,c
        real(dp)::d1,d2,d3
        d1=(x(1)-x(2))*(x(1)-x(3))
        d2=(x(2)-x(1))*(x(2)-x(3))
        d3=(x(3)-x(1))*(x(3)-x(2))
        a=y(1)/d1+y(2)/d2+y(3)/d3
        b=-y(1)*(x(2)+x(3))/d1-y(2)*(x(1)+x(3))/d2-y(3)*(x(1)+x(2))/d3
        c=y(1)*x(2)*x(3)/d1+y(2)*x(1)*x(3)/d2+y(3)*x(1)*x(2)/d3
    end subroutine quadratic_fit
end module boot_profiles
