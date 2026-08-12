! SPDX-License-Identifier: GPL-2.0-only
module clue_dissimilarity
    use clue_kinds, only: dp
    use clue_lsap, only: solve_lsap
    use clue_partition, only: canonicalize_ids, partition_join
    use clue_agreement, only: agreement_rand, contingency_table
    use lpsolve, only: lp_transport, lp_result, LP_MIN, LP_MAX, LP_EQ, LP_OPTIMAL
    implicit none
    private
    public :: dissimilarity_euclidean, dissimilarity_manhattan, dissimilarity_comemberships
    public :: dissimilarity_symdiff, dissimilarity_rand, dissimilarity_gv1
    public :: dissimilarity_ba_a, dissimilarity_ba_c, dissimilarity_ba_d, dissimilarity_ba_e, dissimilarity_vi
    public :: dissimilarity_mallows, dissimilarity_cssd
    public :: hierarchy_dissimilarity_euclidean, hierarchy_dissimilarity_manhattan
    public :: hierarchy_dissimilarity_cophenetic, hierarchy_dissimilarity_gamma
    public :: hierarchy_dissimilarity_chebyshev, hierarchy_dissimilarity_lyapunov
    public :: hierarchy_dissimilarity_spectral
contains
    subroutine conform(mx,my,x,y)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp),allocatable,intent(out)::x(:,:),y(:,:)
        integer::k
        k=max(size(mx,2),size(my,2))
        allocate(x(size(mx,1),k),y(size(my,1),k))
        x=0
        y=0
        if(size(mx,2)>0)x(:,1:size(mx,2))=mx
        if(size(my,2)>0)y(:,1:size(my,2))=my
    end subroutine
    function dissimilarity_euclidean(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        call conform(mx,my,x,y)
        c=matmul(transpose(x),y)
        call solve_lsap(c,p,.true.)
        v=sqrt(sum((x-y(:,p))**2))
    end function
    function dissimilarity_manhattan(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        integer::i,j
        call conform(mx,my,x,y)
        allocate(c(size(x,2),size(y,2)))
        do i=1,size(c,1)
        do j=1,size(c,2)
        c(i,j)=sum(abs(x(:,i)-y(:,j)))
        end do
        end do
        call solve_lsap(c,p)
        v=0
        do i=1,size(p)
        v=v+c(i,p(i))
        end do
    end function
    function dissimilarity_comemberships(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::x(:,:),y(:,:),xx(:,:),xy(:,:),yy(:,:)
        real(dp)::q
        call conform(mx,my,x,y)
        xx=matmul(transpose(x),x)
        xy=matmul(transpose(x),y)
        yy=matmul(transpose(y),y)
        q=sum(xx**2)-2.0_dp*sum(xy**2)+sum(yy**2)
        v=sqrt(max(q,0.0_dp))
    end function
    function dissimilarity_rand(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        v=1.0_dp-agreement_rand(a,b)
    end function
    function dissimilarity_symdiff(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer::n
        n=size(a)
        v=dissimilarity_rand(a,b)*real(n*(n-1),dp)/2.0_dp
    end function
    function dissimilarity_gv1(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::c(:,:),csq(:,:)
        integer,allocatable::p(:)
        integer::kx,ky,k,i
        kx=size(mx,2)
        ky=size(my,2)
        k=max(kx,ky)
        allocate(c(k,k))
        c=0
        if(kx>0.and.ky>0)then
            allocate(csq(kx,ky))
            csq=spread(sum(mx**2,dim=1),2,ky)+spread(sum(my**2,dim=1),1,kx)-2.0_dp*matmul(transpose(mx),my)
            c(1:kx,1:ky)=csq
        end if
        call solve_lsap(c,p)
        v=0
        do i=1,k
        v=v+c(i,p(i))
        end do
        v=sqrt(max(v,0.0_dp))
    end function
    function dissimilarity_ba_a(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        real(dp),allocatable::mx(:,:),my(:,:)
        integer,allocatable::ca(:),cb(:)
        integer::ka,kb
        ca=canonicalize_ids(a)
        cb=canonicalize_ids(b)
        ka=maxval(ca)
        kb=maxval(cb)
        allocate(mx(size(a),ka),my(size(b),kb))
        mx=0
        my=0
        call fillm(ca,mx)
        call fillm(cb,my)
        v=dissimilarity_manhattan(mx,my)/2.0_dp
    contains
        subroutine fillm(ids,m)
        integer,intent(in)::ids(:)
        real(dp),intent(inout)::m(:,:)
        integer::i
        do i=1,size(ids)
        m(i,ids(i))=1
        end do
        end subroutine
    end function
    function dissimilarity_ba_c(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::ca(:),cb(:),j(:)
        ca=canonicalize_ids(a)
        cb=canonicalize_ids(b)
        j=partition_join(ca,cb)
        v=real(maxval(ca)+maxval(cb)-2*maxval(j),dp)
    end function
    function dissimilarity_ba_d(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        v=dissimilarity_rand(a,b)
    end function
    function dissimilarity_ba_e(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp),allocatable::z(:,:),px(:),py(:)
        real(dp)::mi,h,q
        integer::i,j,n
        n=size(a)
        t=contingency_table(a,b)
        allocate(z(size(t,1),size(t,2)))
        z=real(t,dp)/real(n,dp)
        px=sum(z,dim=2)
        py=sum(z,dim=1)
        mi=0
        h=0
        do j=1,size(z,2)
        do i=1,size(z,1)
        if(z(i,j)>0)then
        q=px(i)*py(j)
        if(q>0)mi=mi+z(i,j)*log(z(i,j)/q)
        h=h-z(i,j)*log(z(i,j))
        end if
        end do
        end do
        if(h>0)then
        v=1.0_dp-mi/h
        else
        v=0
        end if
    end function
    function dissimilarity_vi(mx,my,weights) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp),intent(in),optional::weights(:)
        real(dp)::v
        real(dp),allocatable::w(:),px(:),py(:),g(:,:),d(:,:)
        real(dp)::hx,hy
        integer::i,j,n
        n=size(mx,1)
        allocate(w(n))
        if(present(weights))then
        w=weights
        if(size(weights)/=n)w=1
        else
        w=1
        end if
        w=w/sum(w)
        px=matmul(transpose(mx),w)
        py=matmul(transpose(my),w)
        g=matmul(transpose(mx*spread(w,2,size(mx,2))),my)
        allocate(d(size(px),size(py)))
        d=spread(px,2,size(py))*spread(py,1,size(px))
        hx=0
        hy=0
        do i=1,size(px)
        if(px(i)>0)hx=hx-px(i)*log(px(i))
        end do
        do j=1,size(py)
        if(py(j)>0)hy=hy-py(j)*log(py(j))
        end do
        v=hx+hy
        do j=1,size(py)
        do i=1,size(px)
        if(g(i,j)>0.and.d(i,j)>0)v=v-2.0_dp*g(i,j)*log(g(i,j)/d(i,j))
        end do
        end do
    end function
    function dissimilarity_mallows(mx,my,power,alpha,beta,status) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp),intent(in),optional::power,alpha(:),beta(:)
        integer,intent(out),optional::status
        real(dp)::v
        real(dp)::p
        real(dp),allocatable::x(:,:),y(:,:),c(:,:),aa(:),bb(:),flow(:,:)
        integer,allocatable::rs(:),cs(:),none(:)
        type(lp_result)::res
        integer::i,j,k
        p=1
        if(present(power))p=power
        call conform(mx,my,x,y)
        k=size(x,2)
        allocate(c(k,k))
        do i=1,k
        do j=1,k
        c(i,j)=sum(abs(x(:,i)-y(:,j))**p)
        end do
        end do
        allocate(aa(k),bb(k))
        aa=1
        bb=1
        if(present(alpha))aa=alpha
        if(present(beta))bb=beta
        allocate(rs(k),cs(k),flow(k,k),none(0))
        rs=LP_EQ
        cs=LP_EQ
        call lp_transport(c,rs,aa,cs,bb,res,direction=LP_MIN,integer_variables=none,flow=flow)
        if(present(status))status=res%status
        if(res%status==LP_OPTIMAL)then
        v=res%objective**(1.0_dp/p)
        else
        v=huge(1.0_dp)
        end if
    end function
    function dissimilarity_cssd(mx,my,L,alpha,beta,status) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:),L(:,:)
        real(dp),intent(in),optional::alpha(:),beta(:)
        integer,intent(out),optional::status
        real(dp)::v
        real(dp),allocatable::c(:,:),q(:,:),aa(:),bb(:),flow(:,:)
        integer,allocatable::rs(:),cs(:),none(:)
        type(lp_result)::res
        integer::kx,ky
        kx=size(mx,2)
        ky=size(my,2)
        c=matmul(transpose(mx),my)*L
        allocate(aa(kx),bb(ky))
        aa=1
        bb=1
        if(present(alpha))aa=alpha
        if(present(beta))bb=beta
        q=c/(spread(aa,2,ky)+spread(bb,1,kx))
        allocate(rs(kx),cs(ky),flow(kx,ky),none(0))
        rs=LP_EQ
        cs=LP_EQ
        call lp_transport(q,rs,aa,cs,bb,res,direction=LP_MAX,integer_variables=none,flow=flow)
        if(present(status))status=res%status
        if(res%status==LP_OPTIMAL)then
        v=sum(c)-2.0_dp*res%objective
        else
        v=huge(1.0_dp)
        end if
    end function
    function lower_values(a) result(x)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable::x(:)
    integer::i,j,k,n
    n=size(a,1)
    allocate(x(n*(n-1)/2))
    k=0
    do i=2,n
    do j=1,i-1
    k=k+1
    x(k)=a(i,j)
    end do
    end do
    end function
    pure function corr(x,y) result(r)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::r,mx,my,den
    mx=sum(x)/max(1,size(x))
    my=sum(y)/max(1,size(y))
    den=sqrt(sum((x-mx)**2)*sum((y-my)**2))
    if(den>0)then
    r=sum((x-mx)*(y-my))/den
    else
    r=0
    end if
    end function
    function hierarchy_dissimilarity_euclidean(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:)
    x=lower_values(a)
    y=lower_values(b)
    v=sqrt(sum((x-y)**2))
    end function
    function hierarchy_dissimilarity_manhattan(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:)
    x=lower_values(a)
    y=lower_values(b)
    v=sum(abs(x-y))
    end function
    function hierarchy_dissimilarity_cophenetic(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:)
    x=lower_values(a)
    y=lower_values(b)
    v=1.0_dp-corr(x,y)**2
    end function
    function hierarchy_dissimilarity_gamma(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:)
    integer::i,j,cnt,tot
    x=lower_values(a)
    y=lower_values(b)
    cnt=0
    tot=0
    do i=1,size(x)-1
    do j=i+1,size(x)
    tot=tot+1
    if((x(i)-x(j))*(y(i)-y(j))<0)cnt=cnt+1
    end do
    end do
    if(tot>0)then
    v=real(cnt,dp)/tot
    else
    v=0
    end if
    end function
    function hierarchy_dissimilarity_chebyshev(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:)
    x=lower_values(a)
    y=lower_values(b)
    v=maxval(abs(x-y))
    end function
    function hierarchy_dissimilarity_lyapunov(a,b) result(v)
    real(dp),intent(in)::a(:,:),b(:,:)
    real(dp)::v
    real(dp),allocatable::x(:),y(:),q(:)
    x=lower_values(a)
    y=lower_values(b)
    q=x/y
    v=log(maxval(q)/minval(q))
    end function
    function hierarchy_dissimilarity_spectral(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v
        real(dp),allocatable::d(:,:),x(:),y(:)
        integer::it,n
        d=a-b
        n=size(d,1)
        allocate(x(n),y(n))
        x=1.0_dp/sqrt(real(max(1,n),dp))
        v=0
        do it=1,100
        y=matmul(transpose(d),matmul(d,x))
        if(sqrt(sum(y*y))<=tiny(1.0_dp))exit
        x=y/sqrt(sum(y*y))
        end do
        y=matmul(d,x)
        v=sqrt(sum(y*y))
    end function
end module clue_dissimilarity
