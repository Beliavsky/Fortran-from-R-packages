module greybox_optimizer
    use greybox_kinds, only: dp
    implicit none
    private
    public :: pattern_search, numerical_hessian

    abstract interface
        function objective_fn(x) result(f)
            import dp
            real(dp),intent(in)::x(:)
            real(dp)::f
        end function objective_fn
    end interface

contains
    subroutine pattern_search(fun,x,fval,max_iter,tol,lower,upper)
        procedure(objective_fn)::fun
        real(dp),intent(inout)::x(:)
        real(dp),intent(out)::fval
        integer,intent(in),optional::max_iter
        real(dp),intent(in),optional::tol,lower(:),upper(:)
        integer::itmax,it,j
        real(dp)::eps,step(size(x)),candidate(size(x)),fc,best,scale
        logical::improved
        itmax=1000;if(present(max_iter))itmax=max_iter
        eps=1.0e-7_dp;if(present(tol))eps=tol
        do j=1,size(x);step(j)=0.2_dp*max(1.0_dp,abs(x(j)));end do
        if(present(lower))x=max(x,lower);if(present(upper))x=min(x,upper)
        best=fun(x)
        do it=1,itmax
            improved=.false.
            do j=1,size(x)
                candidate=x;candidate(j)=candidate(j)+step(j)
                if(present(upper))candidate(j)=min(candidate(j),upper(j))
                if(present(lower))candidate(j)=max(candidate(j),lower(j))
                fc=fun(candidate)
                if(fc<best)then;x=candidate;best=fc;improved=.true.;cycle;end if
                candidate=x;candidate(j)=candidate(j)-step(j)
                if(present(upper))candidate(j)=min(candidate(j),upper(j))
                if(present(lower))candidate(j)=max(candidate(j),lower(j))
                fc=fun(candidate)
                if(fc<best)then;x=candidate;best=fc;improved=.true.;end if
            end do
            if(.not.improved)step=0.5_dp*step
            scale=max(1.0_dp,maxval(abs(x)))
            if(maxval(step)<=eps*scale)exit
        end do
        fval=best
    end subroutine pattern_search

    subroutine numerical_hessian(fun,x,h)
        procedure(objective_fn)::fun
        real(dp),intent(in)::x(:)
        real(dp),intent(out)::h(:,:)
        real(dp)::xp(size(x)),xm(size(x)),xpp(size(x)),xpm(size(x)),xmp(size(x)),xmm(size(x))
        real(dp)::hi,hj,f0
        integer::i,j
        f0=fun(x);h=0.0_dp
        do i=1,size(x)
            hi=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(i)))
            xp=x;xm=x;xp(i)=xp(i)+hi;xm(i)=xm(i)-hi
            h(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(hi*hi)
            do j=i+1,size(x)
                hj=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x(j)))
                xpp=x;xpm=x;xmp=x;xmm=x
                xpp(i)=xpp(i)+hi;xpp(j)=xpp(j)+hj
                xpm(i)=xpm(i)+hi;xpm(j)=xpm(j)-hj
                xmp(i)=xmp(i)-hi;xmp(j)=xmp(j)+hj
                xmm(i)=xmm(i)-hi;xmm(j)=xmm(j)-hj
                h(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*hi*hj);h(j,i)=h(i,j)
            end do
        end do
    end subroutine numerical_hessian
end module greybox_optimizer
