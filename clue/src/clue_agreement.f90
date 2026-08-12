! SPDX-License-Identifier: GPL-2.0-only
module clue_agreement
    use clue_kinds, only: dp
    use clue_partition, only: canonicalize_ids, membership_from_ids, class_ids_from_membership, partition_join
    use clue_lsap, only: solve_lsap
    implicit none
    private
    public :: agreement_rand, agreement_adjusted_rand, agreement_nmi, agreement_kp
    public :: agreement_angle, agreement_diag, agreement_fm, agreement_jaccard
    public :: agreement_purity, agreement_prediction_strength
    public :: agreement_euclidean, agreement_manhattan
    public :: hierarchy_agreement_euclidean, hierarchy_agreement_manhattan
    public :: hierarchy_agreement_cophenetic, hierarchy_agreement_angle, hierarchy_agreement_gamma
    public :: contingency_table
contains
    function contingency_table(a,b) result(t)
        integer,intent(in)::a(:),b(:)
        integer,allocatable::t(:,:)
        integer,allocatable::ca(:),cb(:)
        integer::i,ka,kb
        if(size(a)/=size(b)) then
        allocate(t(0,0))
        return
        end if
        ca=canonicalize_ids(a)
        cb=canonicalize_ids(b)
        if(size(a)==0) then
        allocate(t(0,0))
        return
        end if
        ka=maxval(ca)
        kb=maxval(cb)
        allocate(t(ka,kb))
        t=0
        do i=1,size(a)
        t(ca(i),cb(i))=t(ca(i),cb(i))+1
        end do
    end function

    pure real(dp) function choose2r(n) result(v)
        integer,intent(in)::n
        v=real(n*(n-1),dp)/2.0_dp
    end function

    function agreement_rand(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp)::sxy,sx,sy
        integer::n
        n=size(a)
        if(n<2 .or. size(b)/=n) then
        v=1.0_dp
        return
        end if
        t=contingency_table(a,b)
        sxy=sum(real(t,dp)**2)
        sx=sum(real(sum(t,dim=2),dp)**2)
        sy=sum(real(sum(t,dim=1),dp)**2)
        v=1.0_dp+(sxy-0.5_dp*(sx+sy))/choose2r(n)
    end function

    function agreement_adjusted_rand(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp)::txy,tx,ty,f,den
        integer::n
        n=size(a)
        if(n<2 .or. size(b)/=n) then
        v=1.0_dp
        return
        end if
        t=contingency_table(a,b)
        txy=sum(real(t,dp)**2)-n
        tx=sum(real(sum(t,dim=2),dp)**2)-n
        ty=sum(real(sum(t,dim=1),dp)**2)-n
        f=tx*ty/real(n*n-n,dp)
        den=0.5_dp*(tx+ty)-f
        if(abs(den)<=epsilon(den)) then
        v=merge(1.0_dp,0.0_dp,abs(txy-f)<=epsilon(den))
        else
        v=(txy-f)/den
        end if
    end function

    function agreement_nmi(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp),allocatable::p(:,:),px(:),py(:)
        real(dp)::mi,ex,ey,q
        integer::i,j,n
        n=size(a)
        if(n==0 .or. size(b)/=n) then
        v=0.0_dp
        return
        end if
        t=contingency_table(a,b)
        allocate(p(size(t,1),size(t,2)))
        p=real(t,dp)/real(n,dp)
        px=sum(p,dim=2)
        py=sum(p,dim=1)
        mi=0.0_dp
        do j=1,size(p,2)
        do i=1,size(p,1)
            q=px(i)*py(j)
            if(p(i,j)>0.0_dp .and. q>0.0_dp) mi=mi+p(i,j)*log(p(i,j)/q)
        end do
        end do
        ex=0.0_dp
        ey=0.0_dp
        do i=1,size(px)
        if(px(i)>0) ex=ex+px(i)*log(px(i))
        end do
        do j=1,size(py)
        if(py(j)>0) ey=ey+py(j)*log(py(j))
        end do
        if(ex*ey<=0.0_dp) then
        v=1.0_dp
        else
        v=mi/sqrt(ex*ey)
        end if
    end function

    function agreement_kp(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp)::axy,ax,ay,num,den,n2
        integer::n
        n=size(a)
        if(n==0 .or. size(b)/=n) then
        v=0
        return
        end if
        t=contingency_table(a,b)
        axy=sum(real(t,dp)**2)
        ax=sum(real(sum(t,dim=2),dp)**2)
        ay=sum(real(sum(t,dim=1),dp)**2)
        n2=real(n,dp)**2
        num=n2*axy-ax*ay
        den=sqrt(max(0.0_dp,ax*(n2-ax)*ay*(n2-ay)))
        if(den<=tiny(1.0_dp)) then
        v=merge(1.0_dp,0.0_dp,abs(num)<=tiny(1.0_dp))
        else
        v=num/den
        end if
    end function

    function agreement_fm(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp)::z,x,y,den
        integer::n
        n=size(a)
        t=contingency_table(a,b)
        z=sum(real(t,dp)**2)-n
        x=sum(real(sum(t,dim=2),dp)**2)-n
        y=sum(real(sum(t,dim=1),dp)**2)-n
        den=sqrt(max(0.0_dp,x*y))
        if(den<=tiny(1.0_dp)) then
        v=merge(1.0_dp,0.0_dp,abs(z)<=tiny(1.0_dp))
        else
        v=z/den
        end if
    end function

    function agreement_jaccard(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        real(dp)::z,den
        integer::n
        n=size(a)
        t=contingency_table(a,b)
        z=sum(real(t,dp)**2)
        den=sum(real(sum(t,dim=2),dp)**2)+sum(real(sum(t,dim=1),dp)**2)-n-z
        if(den<=tiny(1.0_dp)) then
        v=1.0_dp
        else
        v=(z-n)/den
        end if
    end function

    function agreement_purity(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:)
        integer::i,s
        t=contingency_table(a,b)
        s=0
        do i=1,size(t,1)
        if(size(t,2)>0) s=s+maxval(t(i,:))
        end do
        if(size(a)>0) then
        v=real(s,dp)/size(a)
        else
        v=0
        end if
    end function

    function agreement_prediction_strength(a,b) result(v)
        integer,intent(in)::a(:),b(:)
        real(dp)::v
        integer,allocatable::t(:,:),rs(:)
        real(dp)::q
        integer::i
        t=contingency_table(a,b)
        rs=sum(t,dim=2)
        v=huge(1.0_dp)
        do i=1,size(t,1)
            if(rs(i)>1) then
                q=(sum(real(t(i,:),dp)**2)-rs(i))/real(rs(i)*(rs(i)-1),dp)
                v=min(v,q)
            end if
        end do
        if(v>=0.5_dp*huge(1.0_dp)) v=0.0_dp
    end function

    subroutine conform_memberships(mx,my,x,y)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp),allocatable,intent(out)::x(:,:),y(:,:)
        integer::k
        k=max(size(mx,2),size(my,2))
        allocate(x(size(mx,1),k),y(size(my,1),k))
        x=0
        y=0
        if(size(mx,2)>0) x(:,1:size(mx,2))=mx
        if(size(my,2)>0) y(:,1:size(my,2))=my
    end subroutine

    function agreement_angle(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        real(dp)::den
        call conform_memberships(mx,my,x,y)
        c=matmul(transpose(x),y)
        call solve_lsap(c,p,.true.)
        v=0.0_dp
        if(size(p)>0) v=sum(x*y(:,p))
        den=sqrt(sum(x*x)*sum(y*y))
        if(den>0) v=v/den
    end function

    function agreement_diag(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        call conform_memberships(mx,my,x,y)
        c=matmul(transpose(x),y)
        call solve_lsap(c,p,.true.)
        if(size(x,1)>0) then
        v=sum(x*y(:,p))/real(size(x,1),dp)
        else
        v=0
        end if
    end function

    function agreement_euclidean(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v,dmax
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        call conform_memberships(mx,my,x,y)
        c=matmul(transpose(x),y)
        call solve_lsap(c,p,.true.)
        dmax=sqrt(2.0_dp*real(size(x,1),dp))
        if(dmax<=tiny(1.0_dp)) then
        v=1
        else
        v=1.0_dp-sqrt(sum((x-y(:,p))**2))/dmax
        end if
    end function

    function agreement_manhattan(mx,my) result(v)
        real(dp),intent(in)::mx(:,:),my(:,:)
        real(dp)::v,dmax,d
        real(dp),allocatable::x(:,:),y(:,:),c(:,:)
        integer,allocatable::p(:)
        integer::i,j
        call conform_memberships(mx,my,x,y)
        allocate(c(size(x,2),size(y,2)))
        do i=1,size(c,1)
        do j=1,size(c,2)
        c(i,j)=sum(abs(x(:,i)-y(:,j)))
        end do
        end do
        call solve_lsap(c,p)
        d=0
        do i=1,size(p)
        d=d+c(i,p(i))
        end do
        dmax=2.0_dp*real(size(x,1),dp)
        if(dmax<=tiny(1.0_dp)) then
        v=1
        else
        v=1.0_dp-d/dmax
        end if
    end function

    pure function pearson(x,y) result(r)
        real(dp),intent(in)::x(:),y(:)
        real(dp)::r,mx,my,den
        mx=sum(x)/max(1,size(x))
        my=sum(y)/max(1,size(y))
        den=sqrt(sum((x-mx)**2)*sum((y-my)**2))
        if(den>0) then
        r=sum((x-mx)*(y-my))/den
        else
        r=0
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
    function hierarchy_agreement_euclidean(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v
        real(dp),allocatable::x(:),y(:)
        x=lower_values(a)
        y=lower_values(b)
        v=1.0_dp/(1.0_dp+sqrt(sum((x-y)**2)))
    end function
    function hierarchy_agreement_manhattan(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v
        real(dp),allocatable::x(:),y(:)
        x=lower_values(a)
        y=lower_values(b)
        v=1.0_dp/(1.0_dp+sum(abs(x-y)))
    end function
    function hierarchy_agreement_cophenetic(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v
        real(dp),allocatable::x(:),y(:)
        x=lower_values(a)
        y=lower_values(b)
        v=pearson(x,y)
    end function
    function hierarchy_agreement_angle(a,b) result(v)
        real(dp),intent(in)::a(:,:),b(:,:)
        real(dp)::v,den
        real(dp),allocatable::x(:),y(:)
        x=lower_values(a)
        y=lower_values(b)
        den=sqrt(sum(x*x)*sum(y*y))
        if(den>0) then
        v=sum(x*y)/den
        else
        v=0
        end if
    end function
    function hierarchy_agreement_gamma(a,b) result(v)
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
        v=1.0_dp-real(cnt,dp)/tot
        else
        v=1
        end if
    end function
end module clue_agreement
