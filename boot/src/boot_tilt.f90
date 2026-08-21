module boot_tilt
    use boot_kinds, only : dp
    implicit none
    private
    public :: exponential_tilt
contains
    subroutine exponential_tilt(l,theta,t0,strata,p,lambda,info,tol,maxit)
        real(dp),intent(in)::l(:),theta,t0
        integer,intent(in)::strata(:)
        real(dp),intent(out)::p(size(l)),lambda
        integer,intent(out),optional::info
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        real(dp)::target,f,der,lam,tolerance,meanl,mean2,scale,mx
        integer::it,itmax,i,j,ns
        real(dp),allocatable::e(:)
        if(size(strata)/=size(l))error stop "exponential_tilt: strata mismatch"
        tolerance=1.0e-10_dp
        if(present(tol))tolerance=tol
        itmax=100
        if(present(maxit))itmax=maxit
        target=theta-t0
        lam=0.0_dp
        allocate(e(size(l)))
        if(present(info))info=1
        do it=1,itmax
            p=0.0_dp
            do i=1,size(l)
                if(any(strata(1:i-1)==strata(i)))cycle
                ns=count(strata==strata(i))
                scale=real(ns,dp)
                mx=-huge(1.0_dp)
                do j=1,size(l)
                if(strata(j)==strata(i))mx=max(mx,lam*l(j)/scale)
                end do
                do j=1,size(l)
                if(strata(j)==strata(i))e(j)=exp(lam*l(j)/scale-mx)
                end do
                f=sum(e,mask=strata==strata(i))
                do j=1,size(l)
                if(strata(j)==strata(i))p(j)=e(j)/f
                end do
            end do
            meanl=sum(l*p)
            f=meanl-target
            if(abs(f)<=tolerance*max(1.0_dp,abs(target)))then
                lambda=lam
                if(present(info))info=0
                return
            end if
            der=0.0_dp
            do i=1,size(l)
                if(any(strata(1:i-1)==strata(i)))cycle
                ns=count(strata==strata(i))
                scale=real(ns,dp)
                meanl=sum(l*p,mask=strata==strata(i))
                mean2=sum(l*l*p,mask=strata==strata(i))
                der=der+(mean2-meanl*meanl)/scale
            end do
            if(der<=tiny(1.0_dp))exit
            lam=lam-f/der
        end do
        lambda=lam
    end subroutine exponential_tilt
end module boot_tilt
