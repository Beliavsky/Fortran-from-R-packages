module actuar_hachemeister_v02
    use actuar_kinds, only: dp
    implicit none
    private
    public :: hachemeister_result_t, hachemeister_fit

    type :: hachemeister_result_t
        real(dp), allocatable :: collective(:)
        real(dp), allocatable :: individual(:,:)
        real(dp), allocatable :: adjusted(:,:)
        real(dp), allocatable :: credibility(:,:,:)
        real(dp), allocatable :: between_cov(:,:)
        real(dp) :: process_variance = 0.0_dp
        logical :: converged = .false.
        integer :: iterations = 0
    end type hachemeister_result_t
contains
    function hachemeister_fit(ratios,weights,design,tol,maxit) result(res)
        real(dp),intent(in)::ratios(:,:),weights(:,:),design(:,:)
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(hachemeister_result_t)::res
        integer::nc,nt,p,i,j,k,it,mi,eff,nobs
        real(dp)::eps,s2i,denom
        real(dp),allocatable::wmat(:,:,:),cred(:,:,:),a(:,:),csum(:,:),rhs(:,:),sol(:,:), &
            ind(:,:),coll(:),old(:),diff(:),xtwx(:,:),xtwy(:),beta(:),resid(:),tmp(:,:),ident(:,:)
        logical,allocatable::valid(:)
        logical::ok
        nc=size(ratios,1);nt=size(ratios,2);p=size(design,2)
        if(size(weights,1)/=nc .or. size(weights,2)/=nt .or. size(design,1)/=nt .or. p<1) return
        eps=sqrt(epsilon(1.0_dp));if(present(tol)) eps=tol
        mi=100;if(present(maxit)) mi=maxit
        allocate(wmat(p,p,nc),cred(p,p,nc),a(p,p),csum(p,p),rhs(p,1),sol(p,1),ind(p,nc), &
                 coll(p),old(p),diff(p),xtwx(p,p),xtwy(p),beta(p),resid(nt),tmp(p,p),ident(p,p),valid(nc))
        wmat=0.0_dp;cred=0.0_dp;ind=0.0_dp;valid=.false.;res%process_variance=0.0_dp;eff=0
        call identity_matrix(ident)
        do i=1,nc
            xtwx=0.0_dp;xtwy=0.0_dp;nobs=0
            do j=1,nt
                if(weights(i,j)>0.0_dp) then
                    nobs=nobs+1
                    do k=1,p
                        xtwy(k)=xtwy(k)+weights(i,j)*design(j,k)*ratios(i,j)
                        xtwx(k,:)=xtwx(k,:)+weights(i,j)*design(j,k)*design(j,:)
                    end do
                end if
            end do
            if(nobs<p) cycle
            rhs(:,1)=xtwy;call solve_matrix(xtwx,rhs,sol,ok);if(.not.ok) cycle
            beta=sol(:,1);ind(:,i)=beta
            call invert_matrix(xtwx,tmp,ok);if(.not.ok) cycle
            wmat(:,:,i)=tmp;cred(:,:,i)=ident;valid(i)=.true.;eff=eff+1
            resid=0.0_dp
            do j=1,nt
                if(weights(i,j)>0.0_dp) resid(j)=ratios(i,j)-dot_product(design(j,:),beta)
            end do
            if(nobs>p) then
                s2i=0.0_dp
                do j=1,nt
                    if(weights(i,j)>0.0_dp) s2i=s2i+weights(i,j)*resid(j)**2
                end do
                res%process_variance=res%process_variance+s2i/real(nobs-p,dp)
            end if
        end do
        if(eff<2) return
        res%process_variance=res%process_variance/real(eff,dp)
        coll=0.0_dp
        do i=1,nc;if(valid(i)) coll=coll+ind(:,i);end do
        coll=coll/real(eff,dp)
        do it=1,mi
            old=coll;a=0.0_dp
            do i=1,nc
                if(.not.valid(i)) cycle
                diff=ind(:,i)-coll
                a=a+matmul(cred(:,:,i),outer_product(diff,diff))
            end do
            a=0.5_dp*(a+transpose(a))/real(eff-1,dp)
            csum=0.0_dp;rhs=0.0_dp
            do i=1,nc
                if(.not.valid(i)) cycle
                call invert_matrix(a+res%process_variance*wmat(:,:,i),tmp,ok)
                if(.not.ok) cycle
                cred(:,:,i)=matmul(a,tmp)
                csum=csum+cred(:,:,i)
                rhs(:,1)=rhs(:,1)+matmul(cred(:,:,i),ind(:,i))
            end do
            call solve_matrix(csum,rhs,sol,ok);if(.not.ok) exit
            coll=sol(:,1)
            denom=max(1.0_dp,maxval(abs(old)))
            if(maxval(abs(coll-old))/denom<eps) then;res%converged=.true.;exit;end if
        end do
        res%iterations=min(it,mi)
        a=0.0_dp
        do i=1,nc
            if(.not.valid(i)) cycle
            diff=ind(:,i)-coll
            a=a+matmul(cred(:,:,i),outer_product(diff,diff))
        end do
        a=0.5_dp*(a+transpose(a))/real(eff-1,dp)
        do i=1,nc
            if(.not.valid(i)) cycle
            call invert_matrix(a+res%process_variance*wmat(:,:,i),tmp,ok)
            if(ok) cred(:,:,i)=matmul(a,tmp)
        end do
        allocate(res%collective(p),res%individual(p,nc),res%adjusted(p,nc), &
                 res%credibility(p,p,nc),res%between_cov(p,p))
        res%collective=coll;res%individual=ind;res%credibility=cred;res%between_cov=a
        do i=1,nc
            if(valid(i)) then
                res%adjusted(:,i)=coll+matmul(cred(:,:,i),ind(:,i)-coll)
            else
                res%adjusted(:,i)=coll
            end if
        end do
    end function hachemeister_fit

    pure function outer_product(a,b) result(c)
        real(dp),intent(in)::a(:),b(:)
        real(dp),allocatable::c(:,:)
        integer::i,j
        allocate(c(size(a),size(b)))
        do i=1,size(a);do j=1,size(b);c(i,j)=a(i)*b(j);end do;end do
    end function outer_product

    pure subroutine identity_matrix(a)
        real(dp),intent(out)::a(:,:)
        integer::i
        a=0.0_dp;do i=1,min(size(a,1),size(a,2));a(i,i)=1.0_dp;end do
    end subroutine identity_matrix

    pure subroutine invert_matrix(a,ainv,ok)
        real(dp),intent(in)::a(:,:)
        real(dp),allocatable,intent(out)::ainv(:,:)
        logical,intent(out)::ok
        real(dp),allocatable::ident(:,:)
        allocate(ident(size(a,1),size(a,1)));call identity_matrix(ident)
        call solve_matrix(a,ident,ainv,ok)
    end subroutine invert_matrix

    pure subroutine solve_matrix(a,b,x,ok)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp),allocatable,intent(out)::x(:,:)
        logical,intent(out)::ok
        real(dp),allocatable::m(:,:),rhs(:,:),row(:),rrow(:)
        real(dp)::fac
        integer::n,nrhs,i,k,imax
        n=size(a,1);nrhs=size(b,2)
        allocate(m(n,n),rhs(n,nrhs),row(n),rrow(nrhs),x(n,nrhs));m=a;rhs=b;ok=.true.
        do k=1,n
            imax=k
            do i=k+1,n;if(abs(m(i,k))>abs(m(imax,k))) imax=i;end do
            if(abs(m(imax,k))<1.0e-12_dp) then;ok=.false.;x=0.0_dp;return;end if
            if(imax/=k) then
                row=m(k,:);m(k,:)=m(imax,:);m(imax,:)=row
                rrow=rhs(k,:);rhs(k,:)=rhs(imax,:);rhs(imax,:)=rrow
            end if
            do i=k+1,n
                fac=m(i,k)/m(k,k);m(i,k:n)=m(i,k:n)-fac*m(k,k:n);rhs(i,:)=rhs(i,:)-fac*rhs(k,:)
            end do
        end do
        x=0.0_dp
        do i=n,1,-1
            if(i<n) then
                x(i,:)=(rhs(i,:)-matmul(m(i,i+1:n),x(i+1:n,:)))/m(i,i)
            else
                x(i,:)=rhs(i,:)/m(i,i)
            end if
        end do
    end subroutine solve_matrix
end module actuar_hachemeister_v02
