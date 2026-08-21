module boot_saddle
    use boot_kinds, only : dp, pi
    use boot_special, only : normal_cdf, normal_pdf
    use boot_linalg, only : solve_linear, determinant_abs
    implicit none
    private
    public :: multinomial_saddlepoint
contains
    subroutine multinomial_saddlepoint(a,u,mu,pdf,cdf,zeta,strata,lugannani_rice,info)
        real(dp),intent(in)::a(:,:),u(:),mu(:)
        real(dp),intent(out)::pdf,cdf,zeta(size(u))
        integer,intent(in),optional::strata(:)
        logical,intent(in),optional::lugannani_rice
        integer,intent(out),optional::info
        integer::n,d,i,j,k,it,istat,ns
        integer,allocatable::st(:)
        real(dp),allocatable::p0(:),w(:),grad(:),hess(:,:),step(:),meanv(:)
        real(dp)::sw,khat,det2,r,v,eta,mx
        logical::lr
        n=size(a,1)
        d=size(a,2)
        if(size(u)/=d .or. size(mu)/=n)error stop "multinomial_saddlepoint: mismatch"
        allocate(st(n),p0(n),w(n),grad(d),hess(d,d),step(d),meanv(d))
        if(present(strata))then
        st=strata
        else
        st=1
        end if
        do i=1,n
        p0(i)=mu(i)/sum(mu,mask=st==st(i))
        end do
        zeta=0.1_dp
        if(present(info))info=1
        do it=1,100
            grad=-u
            hess=0.0_dp
            do i=1,n
                if(any(st(1:i-1)==st(i)))cycle
                ns=count(st==st(i))
                mx=-huge(1.0_dp)
                do j=1,n
                if(st(j)==st(i))mx=max(mx,dot_product(zeta,a(j,:)))
                end do
                w=0.0_dp
                do j=1,n
                if(st(j)==st(i))w(j)=p0(j)*exp(dot_product(zeta,a(j,:))-mx)
                end do
                sw=sum(w)
                meanv=0.0_dp
                do j=1,n
                if(st(j)==st(i))meanv=meanv+w(j)*a(j,:)/sw
                end do
                grad=grad+real(ns,dp)*meanv
                do j=1,d
                    do k=1,d
                        hess(j,k)=hess(j,k)+real(ns,dp)*(sum(w*a(:,j)*a(:,k))/sw-meanv(j)*meanv(k))
                    end do
                end do
            end do
            call solve_linear(hess,grad,step,istat)
            if(istat/=0)exit
            zeta=zeta-step
            if(maxval(abs(step))<1.0e-10_dp)then
            if(present(info))info=0
            exit
            end if
        end do
        ! recompute cumulant and Hessian at the solution
        khat=-dot_product(zeta,u)
        hess=0.0_dp
        do i=1,n
            if(any(st(1:i-1)==st(i)))cycle
            ns=count(st==st(i))
            mx=-huge(1.0_dp)
            do j=1,n
            if(st(j)==st(i))mx=max(mx,dot_product(zeta,a(j,:)))
            end do
            w=0.0_dp
            do j=1,n
            if(st(j)==st(i))w(j)=p0(j)*exp(dot_product(zeta,a(j,:))-mx)
            end do
            sw=sum(w)
            khat=khat+real(ns,dp)*(log(sw)+mx)
            meanv=0.0_dp
            do j=1,n
            if(st(j)==st(i))meanv=meanv+w(j)*a(j,:)/sw
            end do
            do j=1,d
            do k=1,d
                hess(j,k)=hess(j,k)+real(ns,dp)*(sum(w*a(:,j)*a(:,k))/sw-meanv(j)*meanv(k))
            end do
            end do
        end do
        det2=determinant_abs(hess)
        pdf=exp(khat)/sqrt((2.0_dp*pi)**d*det2)
        cdf=-1.0_dp
        if(d==1)then
            if(abs(zeta(1))<1.0e-10_dp .or. abs(khat)<1.0e-14_dp)then
                cdf=0.5_dp
            else
                r=sign(1.0_dp,zeta(1))*sqrt(max(0.0_dp,-2.0_dp*khat))
                v=zeta(1)*sqrt(det2)
                lr=.false.
                if(present(lugannani_rice))lr=lugannani_rice
                if(lr)then
                    cdf=normal_cdf(r)+normal_pdf(r)*(1.0_dp/r-1.0_dp/v)
                else
                    eta=r+log(v/r)/r
                    cdf=normal_cdf(eta)
                end if
                cdf=max(0.0_dp,min(1.0_dp,cdf))
            end if
        end if
    end subroutine multinomial_saddlepoint
end module boot_saddle
