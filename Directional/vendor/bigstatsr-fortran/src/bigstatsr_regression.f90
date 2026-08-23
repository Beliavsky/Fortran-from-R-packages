! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_regression
    use bigstatsr_kinds, only: dp
    use bigstatsr_fbm, only: fbm_real
    use bigstatsr_utils, only: sigmoid, soft_threshold
    use rspectra_external, only: dgesvd
    implicit none
    private

    type, public :: univ_reg_result
        real(dp), allocatable :: estim(:)
        real(dp), allocatable :: std_err(:)
        real(dp), allocatable :: score(:)
        integer, allocatable :: niter(:)
        logical, allocatable :: converged(:)
    end type univ_reg_result

    type, public :: enet_path_result
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: intercept(:)
        real(dp), allocatable :: loss(:)
        integer, allocatable :: niter(:)
        integer, allocatable :: nactive(:)
        integer :: info = 0
    end type enet_path_result

    type, public :: summary_result
        real(dp), allocatable :: sum_x(:,:)
        real(dp), allocatable :: sum_xx(:,:)
        real(dp), allocatable :: sum_xy(:,:)
    end type summary_result

    public :: big_univ_linreg, big_univ_logreg
    public :: elastic_net_gaussian_path, elastic_net_logistic_path
    public :: predict_enet, big_summaries

    interface
        subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
            import dp
            integer :: n,nrhs,lda,ldb,ipiv(*),info
            real(dp) :: a(lda,*),b(ldb,*)
        end subroutine dgesv
    end interface

contains

    function big_univ_linreg(x,y,covar,rows,cols,thr) result(res)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(in), optional :: covar(:,:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: thr
        type(univ_reg_result) :: res
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: c(:,:),q(:,:),y2(:),col(:),x2(:)
        real(dp) :: tol,beta_num,beta_deno,beta,rss,yss
        integer :: n,j,k
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr)
        if (size(y)/=n) error stop 'big_univ_linreg: y length mismatch'
        if (present(covar)) then
            if (size(covar,1)/=n) error stop 'big_univ_linreg: covariate rows mismatch'
            allocate(c(n,1+size(covar,2)))
            c(:,1)=1.0_dp
            c(:,2:)=covar
        else
            allocate(c(n,1))
            c(:,1)=1.0_dp
        end if
        tol=1.0e-4_dp
        if (present(thr)) tol=thr
        call svd_basis(c,q,tol)
        k=size(q,2)
        allocate(y2(n),col(n),x2(n))
        y2=y
        if (k>0) y2=y-matmul(q,matmul(transpose(q),y))
        yss=dot_product(y2,y2)
        allocate(res%estim(size(cc)),res%std_err(size(cc)),res%score(size(cc)))
        allocate(res%niter(size(cc)),res%converged(size(cc)))
        res%niter=1
        res%converged=.true.
        do j=1,size(cc)
            call x%read_col(cc(j),col,rr)
            x2=col
            if (k>0) x2=col-matmul(q,matmul(transpose(q),col))
            beta_num=dot_product(x2,y2)
            beta_deno=dot_product(x2,x2)
            if (beta_deno<=tiny(1.0_dp) .or. n-k-1<=0) then
                res%estim(j)=0.0_dp
                res%std_err(j)=huge(1.0_dp)
                res%score(j)=0.0_dp
                res%converged(j)=.false.
            else
                beta=beta_num/beta_deno
                rss=max(0.0_dp,yss-beta_num*beta)
                res%estim(j)=beta
                res%std_err(j)=sqrt(rss/(beta_deno*real(n-k-1,dp)))
                if (res%std_err(j)>0.0_dp) then
                    res%score(j)=beta/res%std_err(j)
                else
                    res%score(j)=0.0_dp
                end if
            end if
        end do
    end function big_univ_linreg

    function big_univ_logreg(x,y,covar,rows,cols,tol,maxiter) result(res)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: y(:)
        real(dp), intent(in), optional :: covar(:,:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: tol
        integer, intent(in), optional :: maxiter
        type(univ_reg_result) :: res
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: design(:,:),col(:),coef(:),covmat(:,:)
        real(dp) :: eps
        integer :: mi,j,it,info,p,n
        logical :: conv
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr)
        if (size(y)/=n) error stop 'big_univ_logreg: y length mismatch'
        if (any(y<0.0_dp) .or. any(y>1.0_dp)) error stop 'big_univ_logreg: y must be 0/1'
        p=2
        if (present(covar)) then
            if (size(covar,1)/=n) error stop 'big_univ_logreg: covariate rows mismatch'
            p=p+size(covar,2)
        end if
        allocate(design(n,p),col(n),coef(p),covmat(p,p))
        design(:,2)=1.0_dp
        if (present(covar)) design(:,3:)=covar
        eps=1.0e-8_dp
        if (present(tol)) eps=tol
        mi=20
        if (present(maxiter)) mi=maxiter
        allocate(res%estim(size(cc)),res%std_err(size(cc)),res%score(size(cc)))
        allocate(res%niter(size(cc)),res%converged(size(cc)))
        do j=1,size(cc)
            call x%read_col(cc(j),col,rr)
            design(:,1)=col
            call logistic_irls(design,y,coef,covmat,it,conv,info,eps,mi)
            res%niter(j)=it
            res%converged(j)=conv .and. info==0
            if (res%converged(j)) then
                res%estim(j)=coef(1)
                res%std_err(j)=sqrt(max(0.0_dp,covmat(1,1)))
                if (res%std_err(j)>0.0_dp) then
                    res%score(j)=res%estim(j)/res%std_err(j)
                else
                    res%score(j)=0.0_dp
                end if
            else
                res%estim(j)=0.0_dp
                res%std_err(j)=huge(1.0_dp)
                res%score(j)=0.0_dp
            end if
        end do
    end function big_univ_logreg

    function elastic_net_gaussian_path(x,y,lambda,alpha,center,scale,pf,rows,cols,eps,maxiter) result(res)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: y(:),lambda(:),alpha
        real(dp), intent(in), optional :: center(:),scale(:),pf(:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: eps
        integer, intent(in), optional :: maxiter
        type(enet_path_result) :: res
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: beta(:),r(:),col(:),ctr(:),scl(:),pen(:)
        real(dp) :: l1,l2,z,newb,shift,tol,maxupd
        integer :: l,j,it,n,p,mi
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr); p=size(cc)
        if (size(y)/=n) error stop 'elastic_net_gaussian_path: y mismatch'
        allocate(ctr(p),scl(p),pen(p))
        ctr=0.0_dp; scl=1.0_dp; pen=1.0_dp
        if (present(center)) ctr=center
        if (present(scale)) scl=scale
        if (present(pf)) pen=pf
        if (any(scl<=0.0_dp) .or. any(pen<0.0_dp)) error stop 'elastic_net_gaussian_path: bad scaling/penalty'
        tol=1.0e-5_dp
        if (present(eps)) tol=eps
        mi=1000
        if (present(maxiter)) mi=maxiter
        allocate(res%beta(p,size(lambda)),res%intercept(size(lambda)),res%loss(size(lambda)))
        allocate(res%niter(size(lambda)),res%nactive(size(lambda)))
        allocate(beta(p),r(n),col(n))
        beta=0.0_dp
        r=y-sum(y)/real(n,dp)
        res%beta=0.0_dp
        do l=1,size(lambda)
            l1=lambda(l)*alpha
            l2=lambda(l)-l1
            do it=1,mi
                maxupd=0.0_dp
                do j=1,p
                    call x%read_col(cc(j),col,rr)
                    col=(col-ctr(j))/scl(j)
                    z=dot_product(col,r)/real(n,dp)+beta(j)
                    newb=soft_threshold(z,l1*pen(j),l2*pen(j))
                    shift=newb-beta(j)
                    if (abs(shift)>0.0_dp) then
                        r=r-shift*col
                        beta(j)=newb
                        maxupd=max(maxupd,shift*shift)
                    end if
                end do
                if (maxupd < tol*max(1.0_dp,dot_product(r,r)/real(n,dp))) exit
            end do
            res%beta(:,l)=beta
            res%intercept(l)=sum(y)/real(n,dp)-dot_product(ctr/scl,beta)
            res%loss(l)=dot_product(r,r)
            res%niter(l)=it
            res%nactive(l)=count(abs(beta)>sqrt(epsilon(1.0_dp)))
        end do
    end function elastic_net_gaussian_path

    function elastic_net_logistic_path(x,y,lambda,alpha,center,scale,pf,rows,cols,eps,maxiter) result(res)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: y(:),lambda(:),alpha
        real(dp), intent(in), optional :: center(:),scale(:),pf(:)
        integer, intent(in), optional :: rows(:),cols(:)
        real(dp), intent(in), optional :: eps
        integer, intent(in), optional :: maxiter
        type(enet_path_result) :: res
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: beta(:),eta(:),prob(:),w(:),workr(:),col(:),ctr(:),scl(:),pen(:)
        real(dp) :: l1,l2,u,v,newb,shift,tol,maxupd,si,sumw,intercept,pi0
        integer :: l,j,it,n,p,mi
        call resolve_indices(x,rows,cols,rr,cc)
        n=size(rr); p=size(cc)
        if (size(y)/=n) error stop 'elastic_net_logistic_path: y mismatch'
        if (any(y<0.0_dp) .or. any(y>1.0_dp)) error stop 'elastic_net_logistic_path: y must be 0/1'
        allocate(ctr(p),scl(p),pen(p))
        ctr=0.0_dp; scl=1.0_dp; pen=1.0_dp
        if (present(center)) ctr=center
        if (present(scale)) scl=scale
        if (present(pf)) pen=pf
        tol=1.0e-5_dp
        if (present(eps)) tol=eps
        mi=1000
        if (present(maxiter)) mi=maxiter
        allocate(res%beta(p,size(lambda)),res%intercept(size(lambda)),res%loss(size(lambda)))
        allocate(res%niter(size(lambda)),res%nactive(size(lambda)))
        allocate(beta(p),eta(n),prob(n),w(n),workr(n),col(n))
        beta=0.0_dp
        pi0=min(1.0_dp-1.0e-8_dp,max(1.0e-8_dp,sum(y)/real(n,dp)))
        intercept=log(pi0/(1.0_dp-pi0))
        eta=intercept
        do l=1,size(lambda)
            l1=lambda(l)*alpha
            l2=lambda(l)-l1
            do it=1,mi
                prob=sigmoid(eta)
                w=max(prob*(1.0_dp-prob),1.0e-4_dp)
                workr=(y-prob)/w
                sumw=sum(w)
                si=sum(y-prob)/sumw
                intercept=intercept+si
                eta=eta+si
                workr=workr-si
                maxupd=si*si
                do j=1,p
                    call x%read_col(cc(j),col,rr)
                    col=(col-ctr(j))/scl(j)
                    v=dot_product(w*col,col)/real(n,dp)
                    if (v<=tiny(1.0_dp)) cycle
                    u=dot_product(w*col,workr)/real(n,dp)+v*beta(j)
                    newb=soft_threshold(u,l1*pen(j),l2*pen(j),v)
                    shift=newb-beta(j)
                    if (abs(shift)>0.0_dp) then
                        workr=workr-shift*col
                        eta=eta+shift*col
                        beta(j)=newb
                        maxupd=max(maxupd,shift*shift*v)
                    end if
                end do
                if (maxupd<tol) exit
            end do
            prob=min(1.0_dp-1.0e-15_dp,max(1.0e-15_dp,sigmoid(eta)))
            res%beta(:,l)=beta
            res%intercept(l)=intercept-dot_product(ctr/scl,beta)
            res%loss(l)=-sum(y*log(prob)+(1.0_dp-y)*log(1.0_dp-prob))
            res%niter(l)=it
            res%nactive(l)=count(abs(beta)>sqrt(epsilon(1.0_dp)))
        end do
    end function elastic_net_logistic_path

    function predict_enet(x,beta,intercept,center,scale,rows,cols,logistic) result(pred)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: beta(:),intercept
        real(dp), intent(in), optional :: center(:),scale(:)
        integer, intent(in), optional :: rows(:),cols(:)
        logical, intent(in), optional :: logistic
        real(dp), allocatable :: pred(:)
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        real(dp) :: ctr,scl
        logical :: islog
        integer :: j
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(beta)/=size(cc)) error stop 'predict_enet: beta mismatch'
        allocate(pred(size(rr)),col(size(rr)))
        pred=intercept
        do j=1,size(cc)
            if (abs(beta(j))<=0.0_dp) cycle
            ctr=0.0_dp; scl=1.0_dp
            if (present(center)) ctr=center(j)
            if (present(scale)) scl=scale(j)
            call x%read_col(cc(j),col,rr)
            pred=pred+beta(j)*(col-ctr)/scl
        end do
        islog=.false.
        if (present(logistic)) islog=logistic
        if (islog) pred=sigmoid(pred)
    end function predict_enet

    function big_summaries(x,y,which_set,k,rows,cols) result(res)
        type(fbm_real), intent(in) :: x
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: which_set(:),k
        integer, intent(in), optional :: rows(:),cols(:)
        type(summary_result) :: res
        integer, allocatable :: rr(:),cc(:)
        real(dp), allocatable :: col(:)
        integer :: i,j,g
        call resolve_indices(x,rows,cols,rr,cc)
        if (size(y)/=size(rr) .or. size(which_set)/=size(rr)) error stop 'big_summaries: mismatch'
        allocate(res%sum_x(k,size(cc)),res%sum_xx(k,size(cc)),res%sum_xy(k,size(cc)),col(size(rr)))
        res%sum_x=0.0_dp; res%sum_xx=0.0_dp; res%sum_xy=0.0_dp
        do j=1,size(cc)
            call x%read_col(cc(j),col,rr)
            do i=1,size(rr)
                g=which_set(i)
                if (g<1 .or. g>k) error stop 'big_summaries: group out of bounds'
                res%sum_x(g,j)=res%sum_x(g,j)+col(i)
                res%sum_xx(g,j)=res%sum_xx(g,j)+col(i)*col(i)
                res%sum_xy(g,j)=res%sum_xy(g,j)+col(i)*y(i)
            end do
        end do
    end function big_summaries

    subroutine logistic_irls(a,y,beta,cov,it,converged,info,tol,maxiter)
        real(dp), intent(in) :: a(:,:),y(:),tol
        real(dp), intent(out) :: beta(:),cov(:,:)
        integer, intent(out) :: it,info
        integer, intent(in) :: maxiter
        logical, intent(out) :: converged
        real(dp), allocatable :: eta(:),p(:),w(:),grad(:),h(:,:),rhs(:,:),step(:)
        integer, allocatable :: ipiv(:)
        integer :: q
        q=size(a,2)
        allocate(eta(size(y)),p(size(y)),w(size(y)),grad(q),h(q,q),rhs(q,1),step(q),ipiv(q))
        beta=0.0_dp
        converged=.false.
        info=0
        do it=1,maxiter
            eta=matmul(a,beta)
            p=sigmoid(eta)
            w=max(p*(1.0_dp-p),1.0e-8_dp)
            grad=matmul(transpose(a),y-p)
            h=matmul(transpose(a),spread(w,2,q)*a)
            rhs(:,1)=grad
            call dgesv(q,1,h,q,ipiv,rhs,q,info)
            if (info/=0) return
            step=rhs(:,1)
            beta=beta+step
            if (maxval(abs(step)/(abs(beta)+abs(step)+1.0e-12_dp))<tol) then
                converged=.true.
                exit
            end if
        end do
        eta=matmul(a,beta)
        p=sigmoid(eta)
        w=max(p*(1.0_dp-p),1.0e-8_dp)
        h=matmul(transpose(a),spread(w,2,q)*a)
        call invert_matrix(h,cov,info)
    end subroutine logistic_irls

    subroutine invert_matrix(a,ainv,info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        real(dp), allocatable :: aa(:,:),rhs(:,:)
        integer, allocatable :: ipiv(:)
        integer :: n,i
        integer, intent(out) :: info
        n=size(a,1)
        allocate(aa(n,n),rhs(n,n),ipiv(n))
        aa=a
        rhs=0.0_dp
        do i=1,n
            rhs(i,i)=1.0_dp
        end do
        call dgesv(n,n,aa,n,ipiv,rhs,n,info)
        if (info==0) ainv=rhs
    end subroutine invert_matrix

    subroutine svd_basis(a,q,tol)
        real(dp), intent(in) :: a(:,:),tol
        real(dp), allocatable, intent(out) :: q(:,:)
        real(dp), allocatable :: worka(:,:),sing(:),u(:,:),vt(:,:),work(:)
        real(dp) :: threshold
        integer :: m,n,mn,lwork,info,k
        m=size(a,1)
        n=size(a,2)
        mn=min(m,n)
        allocate(worka(m,n),sing(mn),u(m,mn),vt(1,1))
        worka=a
        lwork=max(1,5*max(m,n))
        allocate(work(lwork))
        call dgesvd('S','N',m,n,worka,m,sing,u,m,vt,1,work,lwork,info)
        if (info/=0) error stop 'svd_basis: LAPACK dgesvd failed'
        threshold=tol*(sqrt(real(m,dp))+sqrt(real(n,dp))-1.0_dp)
        k=count(sing>threshold)
        allocate(q(m,k))
        if (k>0) q=u(:,1:k)
    end subroutine svd_basis

    subroutine resolve_indices(x,rows,cols,rr,cc)
        type(fbm_real), intent(in) :: x
        integer, intent(in), optional :: rows(:),cols(:)
        integer, allocatable, intent(out) :: rr(:),cc(:)
        integer :: i
        if (present(rows)) then
            rr=rows
        else
            allocate(rr(x%nrow)); rr=[(i,i=1,x%nrow)]
        end if
        if (present(cols)) then
            cc=cols
        else
            allocate(cc(x%ncol)); cc=[(i,i=1,x%ncol)]
        end if
    end subroutine resolve_indices

end module bigstatsr_regression
