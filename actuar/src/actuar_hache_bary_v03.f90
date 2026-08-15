module actuar_hache_bary_v03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use actuar_kinds, only: dp
    implicit none
    private
    public :: hache_barycenter_result_t, hachemeister_barycenter_fit

    type :: hache_barycenter_result_t
        real(dp),allocatable::collective(:),individual(:,:),adjusted(:,:)
        real(dp),allocatable::collective_orth(:),individual_orth(:,:),adjusted_orth(:,:)
        real(dp),allocatable::credibility(:,:,:),between_cov(:,:),transition(:,:)
        real(dp)::process_variance=0.0_dp
        logical::converged=.false.
        integer::iterations=0
    end type hache_barycenter_result_t
contains
    function hachemeister_barycenter_fit(ratios,weights,design,iterative,tol,maxit) result(res)
        real(dp),intent(in)::ratios(:,:),weights(:,:),design(:,:)
        logical,intent(in),optional::iterative
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxit
        type(hache_barycenter_result_t)::res
        integer::nc,nt,p,i,j,k,eff,nobs,mi,itmax
        real(dp)::eps,s2,a,zsum,ws,xw,den,prev
        real(dp),allocatable::avew(:),x(:,:),rmat(:,:),ind(:,:),adj(:,:),coll(:),sig2(:), &
            wdiag(:,:),cred(:,:,:),xtwx(:,:),xtwy(:),beta(:),resid(:),rhs(:,:),sol(:,:)
        logical,allocatable::has(:)
        logical::ok,do_iter
        nc=size(ratios,1);nt=size(ratios,2);p=size(design,2)
        if(size(weights,1)/=nc .or. size(weights,2)/=nt .or. size(design,1)/=nt .or. p<1)return
        eps=sqrt(epsilon(1.0_dp));if(present(tol))eps=tol
        mi=100;if(present(maxit))mi=maxit
        do_iter=.false.;if(present(iterative))do_iter=iterative
        allocate(avew(nt),x(nt,p),rmat(p,p),ind(p,nc),adj(p,nc),coll(p),sig2(nc), &
                 wdiag(p,nc),cred(p,p,nc),xtwx(p,p),xtwy(p),beta(p),resid(nt),has(nc))
        avew=0.0_dp
        do j=1,nt
            do i=1,nc
                if(valid_obs(ratios(i,j),weights(i,j)))avew(j)=avew(j)+weights(i,j)
            end do
        end do
        if(sum(avew)<=0.0_dp .or. any(avew<=0.0_dp))return
        avew=avew/sum(avew)
        call weighted_qr(design,avew,x,rmat,ok);if(.not.ok)return
        ind=0.0_dp;sig2=0.0_dp;wdiag=0.0_dp;has=.false.;eff=0
        do i=1,nc
            xtwx=0.0_dp;xtwy=0.0_dp;nobs=0
            do j=1,nt
                if(valid_obs(ratios(i,j),weights(i,j)))then
                    nobs=nobs+1
                    do k=1,p
                        xtwy(k)=xtwy(k)+weights(i,j)*x(j,k)*ratios(i,j)
                    end do
                    call add_outer(xtwx,x(j,:),weights(i,j))
                end if
            end do
            if(nobs>0 .and. sum(weights(i,:),mask=.not.ieee_is_nan(weights(i,:)))>0.0_dp)then
                call solve_vector(xtwx,xtwy,beta,ok)
                if(ok)then
                    ind(:,i)=beta;has(i)=.true.;eff=eff+1;resid=0.0_dp
                    do j=1,nt
                        if(valid_obs(ratios(i,j),weights(i,j))) &
                            resid(j)=ratios(i,j)-dot_product(x(j,:),beta)
                    end do
                    if(nobs>p)then
                        do j=1,nt
                            if(valid_obs(ratios(i,j),weights(i,j))) &
                                sig2(i)=sig2(i)+weights(i,j)*resid(j)**2
                        end do
                        sig2(i)=sig2(i)/real(nobs-p,dp)
                    end if
                end if
            end if
            wdiag(1,i)=sum_valid_weights(weights(i,:))
            do k=2,p
                do j=1,nt
                    if(.not.ieee_is_nan(weights(i,j)) .and. weights(i,j)>0.0_dp) &
                        wdiag(k,i)=wdiag(k,i)+weights(i,j)*x(j,k)**2
                end do
            end do
        end do
        if(eff<2)return
        s2=sum(sig2,mask=has)/real(eff,dp);res%process_variance=s2
        cred=0.0_dp;coll=0.0_dp
        allocate(res%between_cov(p,p));res%between_cov=0.0_dp;itmax=0
        do k=1,p
            ws=sum(wdiag(k,:),mask=has);if(ws<=0.0_dp)cycle
            xw=sum(wdiag(k,:)*ind(k,:),mask=has)/ws
            den=ws*ws-sum(wdiag(k,:)**2,mask=has)
            if(den>0.0_dp)then
                a=ws*(sum(wdiag(k,:)*(ind(k,:)-xw)**2,mask=has)-real(eff-1,dp)*s2)/den
            else
                a=0.0_dp
            end if
            if(do_iter .and. a>0.0_dp .and. weight_range(weights)>sqrt(epsilon(1.0_dp)))then
                do i=1,mi
                    prev=a;zsum=0.0_dp;xw=0.0_dp
                    do j=1,nc
                        if(has(j) .and. wdiag(k,j)>0.0_dp)then
                            zsum=zsum+1.0_dp/(1.0_dp+s2/(wdiag(k,j)*a))
                            xw=xw+ind(k,j)/(1.0_dp+s2/(wdiag(k,j)*a))
                        end if
                    end do
                    if(zsum<=0.0_dp)exit;xw=xw/zsum;a=0.0_dp
                    do j=1,nc
                        if(has(j) .and. wdiag(k,j)>0.0_dp)then
                            zsum=1.0_dp/(1.0_dp+s2/(wdiag(k,j)*prev))
                            a=a+zsum*(ind(k,j)-xw)**2
                        end if
                    end do
                    a=a/real(eff-1,dp);if(abs(a-prev)/max(abs(prev),tiny(1.0_dp))<eps)exit
                end do
                itmax=max(itmax,i)
            end if
            res%between_cov(k,k)=a
            if(a>0.0_dp)then
                zsum=0.0_dp;coll(k)=0.0_dp
                do j=1,nc
                    if(has(j) .and. wdiag(k,j)>0.0_dp)then
                        cred(k,k,j)=1.0_dp/(1.0_dp+s2/(wdiag(k,j)*a))
                        zsum=zsum+cred(k,k,j);coll(k)=coll(k)+cred(k,k,j)*ind(k,j)
                    end if
                end do
                if(zsum>0.0_dp)coll(k)=coll(k)/zsum
            else
                coll(k)=sum(wdiag(k,:)*ind(k,:),mask=has)/ws
            end if
        end do
        do i=1,nc
            if(has(i))then;adj(:,i)=coll+matmul(cred(:,:,i),ind(:,i)-coll)
            else;adj(:,i)=coll;end if
        end do
        allocate(res%transition(p,p));res%transition=rmat
        allocate(res%collective_orth(p),res%individual_orth(p,nc),res%adjusted_orth(p,nc))
        res%collective_orth=coll;res%individual_orth=ind;res%adjusted_orth=adj
        allocate(res%collective(p),res%individual(p,nc),res%adjusted(p,nc),res%credibility(p,p,nc))
        allocate(rhs(p,1))
        rhs(:,1)=coll;call solve_matrix(rmat,rhs,sol,ok);if(ok)res%collective=sol(:,1)
        do i=1,nc
            rhs(:,1)=ind(:,i);call solve_matrix(rmat,rhs,sol,ok);if(ok)res%individual(:,i)=sol(:,1)
            rhs(:,1)=adj(:,i);call solve_matrix(rmat,rhs,sol,ok);if(ok)res%adjusted(:,i)=sol(:,1)
        end do
        res%credibility=cred;res%iterations=itmax;res%converged=.true.
    end function hachemeister_barycenter_fit

    subroutine weighted_qr(design,w,x,rmat,ok)
        real(dp),intent(in)::design(:,:),w(:)
        real(dp),intent(out)::x(:,:),rmat(:,:)
        logical,intent(out)::ok
        real(dp),allocatable::q(:,:),v(:)
        real(dp)::normv
        integer::n,p,i,j
        n=size(design,1);p=size(design,2);allocate(q(n,p),v(n));q=0.0_dp;rmat=0.0_dp;ok=.true.
        do j=1,p
            v=design(:,j)*sqrt(w)
            do i=1,j-1;rmat(i,j)=dot_product(q(:,i),v);v=v-rmat(i,j)*q(:,i);end do
            normv=sqrt(dot_product(v,v));if(normv<1.0e-12_dp)then;ok=.false.;return;end if
            rmat(j,j)=normv;q(:,j)=v/normv
        end do
        do j=1,p;x(:,j)=q(:,j)/sqrt(w);end do
    end subroutine weighted_qr

    pure logical function valid_obs(y,w)
        real(dp),intent(in)::y,w
        valid_obs=.not.ieee_is_nan(y) .and. .not.ieee_is_nan(w) .and. w>0.0_dp
    end function valid_obs

    pure real(dp) function sum_valid_weights(w) result(v)
        real(dp),intent(in)::w(:)
        integer::i
        v=0.0_dp;do i=1,size(w);if(.not.ieee_is_nan(w(i)) .and. w(i)>0.0_dp)v=v+w(i);end do
    end function sum_valid_weights

    pure real(dp) function weight_range(w) result(v)
        real(dp),intent(in)::w(:,:)
        real(dp)::lo,hi
        integer::i,j
        lo=huge(1.0_dp);hi=-huge(1.0_dp)
        do i=1,size(w,1);do j=1,size(w,2)
            if(.not.ieee_is_nan(w(i,j)))then;lo=min(lo,w(i,j));hi=max(hi,w(i,j));end if
        end do;end do
        v=hi-lo
    end function weight_range

    pure subroutine add_outer(a,x,w)
        real(dp),intent(inout)::a(:,:)
        real(dp),intent(in)::x(:),w
        integer::i,j
        do i=1,size(x);do j=1,size(x);a(i,j)=a(i,j)+w*x(i)*x(j);end do;end do
    end subroutine add_outer

    pure subroutine solve_vector(a,b,x,ok)
        real(dp),intent(in)::a(:,:),b(:)
        real(dp),intent(out)::x(:)
        logical,intent(out)::ok
        real(dp),allocatable::bb(:,:),xx(:,:)
        allocate(bb(size(b),1));bb(:,1)=b;call solve_matrix(a,bb,xx,ok);if(ok)x=xx(:,1);if(.not.ok)x=0.0_dp
    end subroutine solve_vector

    pure subroutine solve_matrix(a,b,x,ok)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp),allocatable,intent(out)::x(:,:)
        logical,intent(out)::ok
        real(dp),allocatable::m(:,:),rhs(:,:),row(:),rrow(:)
        real(dp)::fac
        integer::n,nrhs,i,k,imax
        n=size(a,1);nrhs=size(b,2);allocate(m(n,n),rhs(n,nrhs),row(n),rrow(nrhs),x(n,nrhs))
        m=a;rhs=b;ok=.true.
        do k=1,n
            imax=k;do i=k+1,n;if(abs(m(i,k))>abs(m(imax,k)))imax=i;end do
            if(abs(m(imax,k))<1.0e-12_dp)then;ok=.false.;x=0.0_dp;return;end if
            if(imax/=k)then;row=m(k,:);m(k,:)=m(imax,:);m(imax,:)=row
                rrow=rhs(k,:);rhs(k,:)=rhs(imax,:);rhs(imax,:)=rrow;end if
            do i=k+1,n
                fac=m(i,k)/m(k,k);m(i,k:n)=m(i,k:n)-fac*m(k,k:n);rhs(i,:)=rhs(i,:)-fac*rhs(k,:)
            end do
        end do
        x=0.0_dp
        do i=n,1,-1
            if(i<n)then;x(i,:)=(rhs(i,:)-matmul(m(i,i+1:n),x(i+1:n,:)))/m(i,i)
            else;x(i,:)=rhs(i,:)/m(i,i);end if
        end do
    end subroutine solve_matrix
end module actuar_hache_bary_v03
