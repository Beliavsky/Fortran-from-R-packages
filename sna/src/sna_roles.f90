! Role and positional analysis translated from R/sna roles.R.
! Upstream copyright (C) 2004-2024 Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_roles
    use sna_kinds, only : dp, sna_nan, is_missing, sna_eps
    use sna_prep, only : symmetrize
    implicit none
    private
    public :: blockmodel, blockmodel_expand_density, sedist, redist
contains

    function blockmodel(dat,membership,content,diag) result(bm)
        real(dp),intent(in)::dat(:,:,:)
        integer,intent(in)::membership(:)
        character(len=*),intent(in),optional::content
        logical,intent(in),optional::diag
        real(dp),allocatable::bm(:,:,:)
        character(len=24)::ct
        integer::m,n,r,g,p,q,i,j,k,nv,nr,nc
        real(dp),allocatable::vals(:),rs(:),cs(:)
        real(dp)::s
        logical::dg
        m=size(dat,1)
        n=size(dat,2)
        r=maxval(membership)
        ct='density'
        if(present(content))ct=trim(content)
        dg=.false.
        if(present(diag))dg=diag
        allocate(bm(m,r,r))
        bm=sna_nan()
        do g=1,m
        do p=1,r
        do q=1,r
            nv=0
            do i=1,n
            if(membership(i)/=p)cycle
                do j=1,n
                if(membership(j)/=q)cycle
                    if(.not.dg.and.i==j)cycle
                    if(.not.is_missing(dat(g,i,j)))nv=nv+1
                end do
            end do
            if(nv==0)cycle
            allocate(vals(nv))
            k=0
            do i=1,n
            if(membership(i)/=p)cycle
                do j=1,n
                if(membership(j)/=q)cycle
                    if(.not.dg.and.i==j)cycle
                    if(.not.is_missing(dat(g,i,j)))then
                    k=k+1
                    vals(k)=dat(g,i,j)
                    end if
                end do
            end do
            select case(ct)
            case('density')
            bm(g,p,q)=sum(vals)/real(nv,dp)
            case('sum')
            bm(g,p,q)=sum(vals)
            case('min')
            bm(g,p,q)=minval(vals)
            case('max')
            bm(g,p,q)=maxval(vals)
            case('median')
            bm(g,p,q)=median(vals)
            case('meanrowsum')
                nr=count(membership==p)
                allocate(rs(nr))
                rs=0
                s=0
                k=0
                do i=1,n
                    if(membership(i)/=p)cycle
                    k=k+1
                    do j=1,n
                        if(membership(j)==q.and.(dg.or.i/=j).and..not.is_missing(dat(g,i,j)))rs(k)=rs(k)+dat(g,i,j)
                    end do
                end do
                bm(g,p,q)=sum(rs)/real(nr,dp)
                deallocate(rs)
            case('meancolsum')
                nc=count(membership==q)
                allocate(cs(nc))
                cs=0
                k=0
                do j=1,n
                    if(membership(j)/=q)cycle
                    k=k+1
                    do i=1,n
                        if(membership(i)==p.and.(dg.or.i/=j).and..not.is_missing(dat(g,i,j)))cs(k)=cs(k)+dat(g,i,j)
                    end do
                end do
                bm(g,p,q)=sum(cs)/real(nc,dp)
                deallocate(cs)
            case default
                bm(g,p,q)=sum(vals)/real(nv,dp)
            end select
            deallocate(vals)
        end do
        end do
        end do
    end function blockmodel

    function blockmodel_expand_density(bm,ev,mode,diag) result(expanded)
        real(dp),intent(in)::bm(:,:,:)
        integer,intent(in)::ev(:)
        character(len=*),intent(in),optional::mode
        logical,intent(in),optional::diag
        real(dp),allocatable::expanded(:,:,:)
        integer::m,en,g,p,q,i,j,sp,ep,sq,eq
        real(dp)::u,pr
        logical::dg,undir
        character(len=16)::md
        m=size(bm,1)
        en=sum(ev)
        allocate(expanded(m,en,en))
        expanded=0.0_dp
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        do g=1,m
            sp=1
            do p=1,size(ev)
                ep=sp+ev(p)-1
                sq=1
                do q=1,size(ev)
                    eq=sq+ev(q)-1
                    pr=bm(g,p,q)
                    if(is_missing(pr))pr=0.0_dp
                    do i=sp,ep
                    do j=sq,eq
                        if(.not.dg.and.i==j)cycle
                        if(undir.and.j>i)cycle
                        call random_number(u)
                        if(u<max(0.0_dp,min(1.0_dp,pr)))expanded(g,i,j)=1.0_dp
                        if(undir)expanded(g,j,i)=expanded(g,i,j)
                    end do
                    end do
                    sq=eq+1
                end do
                sp=ep+1
            end do
        end do
    end function blockmodel_expand_density

    function sedist(dat,method,joint_analysis,mode,diag,code_diss) result(o)
        real(dp),intent(in)::dat(:,:,:)
        character(len=*),intent(in),optional::method,mode
        logical,intent(in),optional::joint_analysis,diag,code_diss
        real(dp),allocatable::o(:,:,:)
        real(dp),allocatable::v(:,:),one(:,:,:)
        character(len=24)::meth,md
        logical::joint,dg,diss
        integer::m,n,k,i,j,rows,pos,g
        m=size(dat,1)
        n=size(dat,2)
        meth='hamming'
        if(present(method))meth=trim(method)
        md='digraph'
        if(present(mode))md=trim(mode)
        joint=.false.
        if(present(joint_analysis))joint=joint_analysis
        dg=.false.
        if(present(diag))dg=diag
        diss=.false.
        if(present(code_diss))diss=code_diss
        if(joint)then
            allocate(o(1,n,n))
            rows=2*m*n
            allocate(v(rows,n))
            pos=0
            do g=1,m
                do k=1,n
                    pos=pos+1
                    v(pos,:)=dat(g,:,k) ! incoming profiles
                end do
                do k=1,n
                    pos=pos+1
                    v(pos,:)=dat(g,k,:) ! outgoing profiles, columns are positions
                end do
            end do
            call profile_distance(v,meth,diss,o(1,:,:))
        else
            allocate(o(m,n,n))
            do g=1,m
                rows=2*n
                allocate(v(rows,n))
                pos=0
                do k=1,n
                pos=pos+1
                v(pos,:)=dat(g,:,k)
                end do
                do k=1,n
                pos=pos+1
                v(pos,:)=dat(g,k,:)
                end do
                if(.not.dg)then
                    ! R masks diagonals before profile construction.
                    do i=1,n
                        v(i,i)=sna_nan()
                        v(n+i,i)=sna_nan()
                    end do
                end if
                call profile_distance(v,meth,diss,o(g,:,:))
                deallocate(v)
            end do
        end if
    end function sedist

    subroutine profile_distance(v,method,diss,out)
        real(dp),intent(in)::v(:,:)
        character(len=*),intent(in)::method
        logical,intent(in)::diss
        real(dp),intent(out)::out(:,:)
        integer::i,j,k,c,con,dis
        real(dp)::sx,sy,sxx,syy,sxy,z
        logical::different
        do i=1,size(v,2)
        do j=1,size(v,2)
            select case(trim(method))
            case('euclidean')
                z=0
                do k=1,size(v,1)
                if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))z=z+(v(k,i)-v(k,j))**2
                end do
                out(i,j)=sqrt(z)
            case('hamming')
                z=0
                do k=1,size(v,1)
                if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))z=z+abs(v(k,i)-v(k,j))
                end do
                out(i,j)=z
            case('gamma')
                con=0
                dis=0
                do k=1,size(v,1)
                if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))then
                if(v(k,i)==v(k,j))then
                con=con+1
                else
                dis=dis+1
                end if
                end if
                end do
                if(con+dis>0)then
                out(i,j)=real(con-dis,dp)/real(con+dis,dp)
                else
                out(i,j)=0
                end if
                if(diss)out(i,j)=-out(i,j)
            case('exact')
                different=.false.
                do k=1,size(v,1)
                if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))then
                if(v(k,i)/=v(k,j))different=.true.
                end if
                end do
                out(i,j)=merge(1.0_dp,0.0_dp,different)
            case('correlation')
                sx=0
                sy=0
                c=0
                do k=1,size(v,1)
                if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))then
                sx=sx+v(k,i)
                sy=sy+v(k,j)
                c=c+1
                end if
                end do
                if(c<2)then
                out(i,j)=0
                else
                    sx=sx/real(c,dp)
                    sy=sy/real(c,dp)
                    sxx=0
                    syy=0
                    sxy=0
                    do k=1,size(v,1)
                    if(.not.is_missing(v(k,i)).and..not.is_missing(v(k,j)))then
                    sxx=sxx+(v(k,i)-sx)**2
                    syy=syy+(v(k,j)-sy)**2
                    sxy=sxy+(v(k,i)-sx)*(v(k,j)-sy)
                    end if
                    end do
                    if(sxx<=sna_eps.or.syy<=sna_eps)then
                    out(i,j)=0
                    else
                    out(i,j)=sxy/sqrt(sxx*syy)
                    end if
                end if
                if(diss)out(i,j)=-out(i,j)
            case default
                out(i,j)=0
            end select
        end do
        end do
    end subroutine profile_distance

    function redist(dat,code_diss,seed_partition,mode,diag) result(eq)
        real(dp),intent(in)::dat(:,:,:)
        logical,intent(in),optional::code_diss,diag
        integer,intent(in),optional::seed_partition(:)
        character(len=*),intent(in),optional::mode
        real(dp),allocatable::eq(:,:)
        integer,allocatable::part1(:),part2(:),history(:,:),tmp_hist(:,:)
        integer,allocatable::cat(:,:),signature(:,:,:)
        integer::m,n,i,j,g,bits,r,iter,maxiter,a,b,p,code,nh
        logical::diss,dg,changed,same,undir
        character(len=16)::md
        real(dp)::mx,mn
        m=size(dat,1)
        n=size(dat,2)
        bits=2*m
        r=2**bits-1
        diss=.true.
        if(present(code_diss))diss=code_diss
        dg=.false.
        if(present(diag))dg=diag
        md='digraph'
        if(present(mode))md=trim(mode)
        undir=trim(md)=='graph'
        allocate(cat(n,n))
        cat=0
        do i=1,n
        do j=1,n
            if(.not.dg.and.i==j)cycle
            code=0
            do g=1,m
                if(dat(g,i,j)>0.0_dp)code=ibset(code,2*g-2)
                if(dat(g,j,i)>0.0_dp)code=ibset(code,2*g-1)
            end do
            cat(i,j)=code
        end do
        end do
        allocate(part1(n),part2(n))
        if(present(seed_partition))then
        part1=seed_partition
        else
        part1=1
        end if
        maxiter=n+1
        allocate(history(maxiter,n))
        history=0
        nh=0
        do
            nh=nh+1
            history(nh,:)=part1
            allocate(signature(r,n,n))
            signature=0
            do i=1,n
            do j=1,n
                code=cat(i,j)
                if(code>0)signature(code,i,part1(j))=1
            end do
            end do
            part2=[(i,i=1,n)]
            changed=.false.
            do i=2,n
                do j=1,i-1
                    if(part1(i)/=part1(j))cycle
                    same=all(signature(:,i,:)==signature(:,j,:))
                    if(same)then
                    part2(i)=part2(j)
                    else
                    changed=.true.
                    end if
                end do
            end do
            deallocate(signature)
            part1=part2
            if(.not.changed.or.nh>=maxiter)exit
        end do
        allocate(eq(n,n))
        eq=0
        do i=1,n
        do j=1,n
            do iter=1,nh
                if(history(iter,i)==history(iter,j))eq(i,j)=real(iter,dp)
            end do
        end do
        end do
        if(diss)then
            mx=maxval(eq)
            mn=minval(eq)
            if(mx>mn)eq=(mx-eq)/(mx-mn)
            if(mx<=mn)eq=0
        end if
    end function redist

    function median(x) result(v)
        real(dp),intent(in)::x(:)
        real(dp)::v
        real(dp),allocatable::y(:)
        integer::i,j,n
        real(dp)::t
        n=size(x)
        allocate(y(n))
        y=x
        do i=2,n
        t=y(i)
        j=i-1
        do while(j>=1)
        if(y(j)<=t)exit
        y(j+1)=y(j)
        j=j-1
        end do
        y(j+1)=t
        end do
        if(mod(n,2)==1)then
        v=y((n+1)/2)
        else
        v=0.5_dp*(y(n/2)+y(n/2+1))
        end if
    end function median
end module sna_roles
