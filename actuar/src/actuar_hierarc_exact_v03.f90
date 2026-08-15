module actuar_hierarc_exact_v03
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use actuar_kinds, only: dp
    implicit none
    private
    public :: hierarc_exact_level_t, hierarc_exact_result_t, hierarc_exact_fit
    public :: HIERARC_BUHLMANN_GISLER, HIERARC_OHLSSON, HIERARC_ITERATIVE
    integer,parameter::HIERARC_BUHLMANN_GISLER=1,HIERARC_OHLSSON=2,HIERARC_ITERATIVE=3

    type :: hierarc_exact_level_t
        real(dp),allocatable::weight(:),mean(:),credibility(:),premium(:)
        integer,allocatable::parent(:)
    end type hierarc_exact_level_t
    type :: hierarc_exact_result_t
        integer::nlevels=0,method=0,iterations=0
        real(dp),allocatable::variance(:)
        type(hierarc_exact_level_t),allocatable::level(:)
        logical::converged=.false.
    end type hierarc_exact_result_t
contains
    function hierarc_exact_fit(ratios,weights,classification,method,tol,maxit) result(res)
        real(dp),intent(in)::ratios(:,:),weights(:,:)
        integer,intent(in)::classification(:,:)
        integer,intent(in),optional::method,maxit
        real(dp),intent(in),optional::tol
        type(hierarc_exact_result_t)::res
        integer::nr,np,l,meth,mi,k,i,j,g,p,c,nobs,neff,it,idx
        integer,allocatable::gid(:,:),gcount(:),effnodes(:),effchild(:)
        real(dp),allocatable::roww(:),rowm(:),denom(:),b(:),oldb(:),bi(:),ci(:),childw(:)
        real(dp)::s2,between,num,den,eps,diff,bw,wuse
        nr=size(ratios,1);np=size(ratios,2);l=size(classification,2)
        if(size(weights,1)/=nr .or. size(weights,2)/=np .or. size(classification,1)/=nr .or. l<1)return
        meth=HIERARC_BUHLMANN_GISLER;if(present(method))meth=method
        mi=100;if(present(maxit))mi=maxit;eps=sqrt(epsilon(1.0_dp));if(present(tol))eps=tol
        allocate(gid(l,nr),gcount(l));call build_group_ids(classification,gid,gcount)
        allocate(roww(nr),rowm(nr));roww=0.0_dp;rowm=0.0_dp;nobs=0
        do i=1,nr
            do j=1,np
                if(valid_obs(ratios(i,j),weights(i,j)))then
                    roww(i)=roww(i)+weights(i,j);rowm(i)=rowm(i)+weights(i,j)*ratios(i,j);nobs=nobs+1
                end if
            end do
            if(roww(i)>0.0_dp)rowm(i)=rowm(i)/roww(i)
        end do
        allocate(effnodes(l));effnodes=0
        do k=1,l
            do g=1,gcount(k)
                if(any([(roww(i)>0.0_dp .and. gid(k,i)==g,i=1,nr)]))effnodes(k)=effnodes(k)+1
            end do
        end do
        allocate(denom(l+1));denom(1)=real(effnodes(1)-1,dp)
        do k=2,l;denom(k)=real(effnodes(k)-effnodes(k-1),dp);end do
        denom(l+1)=real(nobs-effnodes(l),dp)
        if(any(denom<=0.0_dp))return
        s2=0.0_dp
        do i=1,nr;do j=1,np
            if(valid_obs(ratios(i,j),weights(i,j)))s2=s2+weights(i,j)*(ratios(i,j)-rowm(i))**2
        end do;end do
        s2=s2/denom(l+1)
        res%nlevels=l;res%method=meth;allocate(res%variance(l+1),res%level(0:l));allocate(b(l+1));b=0.0_dp;b(l+1)=s2
        allocate(res%level(l)%weight(gcount(l)),res%level(l)%mean(gcount(l)))
        res%level(l)%weight=0.0_dp;res%level(l)%mean=0.0_dp
        do i=1,nr
            g=gid(l,i);res%level(l)%weight(g)=res%level(l)%weight(g)+roww(i)
            res%level(l)%mean(g)=res%level(l)%mean(g)+roww(i)*rowm(i)
        end do
        do g=1,gcount(l)
            if(res%level(l)%weight(g)>0.0_dp)res%level(l)%mean(g)=res%level(l)%mean(g)/res%level(l)%weight(g)
        end do
        do k=l,1,-1
            c=gcount(k)
            if(k>1)then;p=gcount(k-1);else;p=1;end if
            allocate(res%level(k)%parent(c));call node_parents(k,gid,gcount,res%level(k)%parent)
            if(.not.allocated(res%level(k-1)%weight))allocate(res%level(k-1)%weight(p),res%level(k-1)%mean(p))
            res%level(k-1)%weight=0.0_dp;res%level(k-1)%mean=0.0_dp
            do i=1,c
                g=res%level(k)%parent(i);res%level(k-1)%weight(g)=res%level(k-1)%weight(g)+res%level(k)%weight(i)
                res%level(k-1)%mean(g)=res%level(k-1)%mean(g)+res%level(k)%weight(i)*res%level(k)%mean(i)
            end do
            do g=1,p
                if(res%level(k-1)%weight(g)>0.0_dp) &
                    res%level(k-1)%mean(g)=res%level(k-1)%mean(g)/res%level(k-1)%weight(g)
            end do
            between=next_nonzero_variance(b,k+1,l+1)
            allocate(bi(p),ci(p),effchild(p));bi=0.0_dp;ci=0.0_dp;effchild=0
            do i=1,c
                g=res%level(k)%parent(i)
                if(res%level(k)%weight(i)>0.0_dp)effchild(g)=effchild(g)+1
                bi(g)=bi(g)+res%level(k)%weight(i)*(res%level(k)%mean(i)-res%level(k-1)%mean(g))**2
                ci(g)=ci(g)+res%level(k)%weight(i)**2
            end do
            do g=1,p
                bi(g)=bi(g)-real(max(0,effchild(g)-1),dp)*between
                if(res%level(k-1)%weight(g)>0.0_dp)then
                    ci(g)=res%level(k-1)%weight(g)-ci(g)/res%level(k-1)%weight(g)
                else;ci(g)=0.0_dp;end if
            end do
            if(meth==HIERARC_BUHLMANN_GISLER)then
                num=0.0_dp;idx=0
                do g=1,p
                    if(ci(g)/=0.0_dp)then;num=num+max(bi(g)/ci(g),0.0_dp);idx=idx+1;end if
                end do
                if(idx>0)b(k)=num/real(idx,dp)
            else
                den=sum(ci);if(den/=0.0_dp)b(k)=sum(bi)/den
            end if
            allocate(res%level(k)%credibility(c));res%level(k)%credibility=0.0_dp
            if(b(k)/=0.0_dp)then
                do i=1,c
                    if(res%level(k)%weight(i)>0.0_dp) &
                        res%level(k)%credibility(i)=1.0_dp/(1.0_dp+between/(b(k)*res%level(k)%weight(i)))
                end do
                call aggregate_with_cred(k,res)
            end if
            deallocate(bi,ci,effchild)
        end do
        if(meth==HIERARC_ITERATIVE)then
            b=max(b,0.0_dp)
            if(any(b(:l)>0.0_dp))then
                allocate(oldb(l+1));oldb(l+1)=b(l+1)
                do it=1,mi
                    oldb(:l)=b(:l);diff=0.0_dp
                    do k=l,1,-1
                        c=gcount(k)
                        if(k>1)then;p=gcount(k-1);else;p=1;end if
                        bw=next_nonzero_variance(b,k+1,l+1)
                        call iterative_aggregate(k,b(k),bw,res)
                        if(oldb(k)>0.0_dp)then
                            b(k)=0.0_dp
                            do i=1,c
                                g=res%level(k)%parent(i)
                                wuse=effective_child_weight(res%level(k)%credibility(i),res%level(k)%weight(i))
                                b(k)=b(k)+wuse*(res%level(k)%mean(i)-res%level(k-1)%mean(g))**2
                            end do
                            b(k)=b(k)/denom(k);if(b(k)<=eps*eps)b(k)=0.0_dp
                        end if
                        call iterative_aggregate(k,b(k),bw,res)
                    end do
                    do k=1,l
                        if(b(k)>0.0_dp .and. oldb(k)>0.0_dp) &
                            diff=max(diff,abs(b(k)-oldb(k))/oldb(k))
                    end do
                    if(diff<eps)exit
                end do
                res%iterations=min(it,mi);res%converged=(it<=mi)
            else
                res%iterations=0;res%converged=.true.
            end if
        else
            res%iterations=0;res%converged=.true.
        end if
        res%variance=b
        allocate(res%level(0)%premium(1));res%level(0)%premium=res%level(0)%mean
        do k=1,l
            c=gcount(k);allocate(res%level(k)%premium(c))
            do i=1,c
                g=res%level(k)%parent(i)
                res%level(k)%premium(i)=res%level(k-1)%premium(g)+res%level(k)%credibility(i)* &
                    (res%level(k)%mean(i)-res%level(k-1)%premium(g))
            end do
        end do
    end function hierarc_exact_fit

    subroutine aggregate_with_cred(k,res)
        integer,intent(in)::k
        type(hierarc_exact_result_t),intent(inout)::res
        integer::i,g
        res%level(k-1)%weight=0.0_dp;res%level(k-1)%mean=0.0_dp
        do i=1,size(res%level(k)%weight)
            g=res%level(k)%parent(i);res%level(k-1)%weight(g)=res%level(k-1)%weight(g)+res%level(k)%credibility(i)
            res%level(k-1)%mean(g)=res%level(k-1)%mean(g)+res%level(k)%credibility(i)*res%level(k)%mean(i)
        end do
        do g=1,size(res%level(k-1)%weight)
            if(res%level(k-1)%weight(g)>0.0_dp) &
                res%level(k-1)%mean(g)=res%level(k-1)%mean(g)/res%level(k-1)%weight(g)
        end do
    end subroutine aggregate_with_cred

    subroutine iterative_aggregate(k,bk,bw,res)
        integer,intent(in)::k
        real(dp),intent(in)::bk,bw
        type(hierarc_exact_result_t),intent(inout)::res
        integer::i,g
        real(dp)::z,wuse
        res%level(k-1)%weight=0.0_dp;res%level(k-1)%mean=0.0_dp
        do i=1,size(res%level(k)%weight)
            if(bk>0.0_dp .and. res%level(k)%weight(i)>0.0_dp)then
                z=1.0_dp/(1.0_dp+bw/(bk*res%level(k)%weight(i)))
            else;z=0.0_dp;end if
            res%level(k)%credibility(i)=z;wuse=effective_child_weight(z,res%level(k)%weight(i))
            g=res%level(k)%parent(i);res%level(k-1)%weight(g)=res%level(k-1)%weight(g)+wuse
            res%level(k-1)%mean(g)=res%level(k-1)%mean(g)+wuse*res%level(k)%mean(i)
        end do
        do g=1,size(res%level(k-1)%weight)
            if(res%level(k-1)%weight(g)>0.0_dp) &
                res%level(k-1)%mean(g)=res%level(k-1)%mean(g)/res%level(k-1)%weight(g)
        end do
    end subroutine iterative_aggregate

    pure real(dp) function effective_child_weight(z,w) result(v)
        real(dp),intent(in)::z,w
        if(z/=0.0_dp)then;v=z;else;v=w;end if
    end function effective_child_weight

    pure real(dp) function next_nonzero_variance(b,first,last) result(v)
        real(dp),intent(in)::b(:)
        integer,intent(in)::first,last
        integer::i
        v=0.0_dp
        do i=first,last;if(b(i)/=0.0_dp)then;v=b(i);return;end if;end do
    end function next_nonzero_variance

    subroutine build_group_ids(classification,gid,gcount)
        integer,intent(in)::classification(:,:)
        integer,intent(out)::gid(:,:),gcount(:)
        integer::nr,l,k,i,j,g
        nr=size(classification,1);l=size(classification,2);gid=0;gcount=0
        do k=1,l
            do i=1,nr
                g=0
                do j=1,i-1
                    if(all(classification(i,1:k)==classification(j,1:k)))then;g=gid(k,j);exit;end if
                end do
                if(g==0)then;gcount(k)=gcount(k)+1;g=gcount(k);end if
                gid(k,i)=g
            end do
        end do
    end subroutine build_group_ids

    subroutine node_parents(k,gid,gcount,parent)
        integer,intent(in)::k,gid(:,:),gcount(:)
        integer,intent(out)::parent(:)
        integer::i,c
        parent=0
        do i=1,size(gid,2)
            c=gid(k,i)
            if(parent(c)==0)then
                if(k==1)then;parent(c)=1;else;parent(c)=gid(k-1,i);end if
            end if
        end do
    end subroutine node_parents

    pure logical function valid_obs(y,w)
        real(dp),intent(in)::y,w
        valid_obs=.not.ieee_is_nan(y) .and. .not.ieee_is_nan(w) .and. w>0.0_dp
    end function valid_obs
end module actuar_hierarc_exact_v03
