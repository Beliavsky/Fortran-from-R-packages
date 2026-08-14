! Core connectivity, census, and graph-level measures translated from R/sna.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_graph
    use sna_kinds, only : dp, sna_inf, sna_nan, is_missing, sna_eps
    use sna_types, only : geodist_result, component_result, path_census_result
    use sna_prep, only : symmetrize, nties
    implicit none
    private

    public :: edge_present, geodist, reachability, component_dist, components
    public :: component_size_byvertex, component_largest_mask, is_connected
    public :: isolates, is_isolate, cutpoints, bicomponent_dist, kcores
    public :: maxflow, neighborhood, simmelian, structure_statistics
    public :: degree, dyad_census, gden, connectedness, efficiency, grecip
    public :: mutuality, triad_classify, triad_census, gtrans, hierarchy, lubness
    public :: kpath_census, kcycle_census, clique_census

contains

    elemental logical function edge_present(x) result(ans)
        real(dp), intent(in) :: x
        ans = (.not.is_missing(x)) .and. x /= 0.0_dp
    end function edge_present

    function geodist(g, ignore_eval, count_paths) result(res)
        real(dp), intent(in) :: g(:,:)
        logical, intent(in), optional :: ignore_eval, count_paths
        type(geodist_result) :: res
        logical :: ign, cp
        integer :: n, s

        n=size(g,1)
        if(size(g,2)/=n) error stop 'geodist: square matrix required'
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        cp=.true.
        if(present(count_paths))cp=count_paths
        allocate(res%distance(n,n),res%counts(n,n))
        res%distance=sna_inf()
        res%counts=0.0_dp
        do s=1,n
            if(ign)then
                call sssp_unweighted(g,s,res%distance(s,:),res%counts(s,:))
            else
                call sssp_weighted(g,s,res%distance(s,:),res%counts(s,:))
            end if
        end do
        if(.not.cp)res%counts=0.0_dp
    end function geodist

    subroutine sssp_unweighted(g,s,d,sigma)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::s
        real(dp),intent(out)::d(:),sigma(:)
        integer,allocatable::q(:)
        integer::head,tail,v,w,n
        n=size(g,1)
        allocate(q(n))
        d=sna_inf()
        sigma=0.0_dp
        d(s)=0.0_dp
        sigma(s)=1.0_dp
        head=1
        tail=1
        q(1)=s
        do while(head<=tail)
            v=q(head)
            head=head+1
            do w=1,n
                if(w==v .or. .not.edge_present(g(v,w)))cycle
                if(d(w)==sna_inf())then
                    d(w)=d(v)+1.0_dp
                    tail=tail+1
                    q(tail)=w
                end if
                if(abs(d(w)-(d(v)+1.0_dp))<=sna_eps)sigma(w)=sigma(w)+sigma(v)
            end do
        end do
    end subroutine sssp_unweighted

    subroutine sssp_weighted(g,s,d,sigma)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::s
        real(dp),intent(out)::d(:),sigma(:)
        logical,allocatable::done(:)
        integer::n,it,v,w,i
        real(dp)::best,alt,wt
        n=size(g,1)
        allocate(done(n))
        done=.false.
        d=sna_inf()
        sigma=0.0_dp
        d(s)=0.0_dp
        sigma(s)=1.0_dp
        do it=1,n
            v=0
            best=sna_inf()
            do i=1,n
                if(.not.done(i).and.d(i)<best)then
                best=d(i)
                v=i
                end if
            end do
            if(v==0)exit
            done(v)=.true.
            do w=1,n
                if(w==v.or..not.edge_present(g(v,w)))cycle
                wt=g(v,w)
                if(wt<0.0_dp)error stop 'geodist: negative edge values not supported'
                alt=d(v)+wt
                if(alt<d(w)-sna_eps)then
                    d(w)=alt
                    sigma(w)=sigma(v)
                else if(abs(alt-d(w))<=sna_eps)then
                    sigma(w)=sigma(w)+sigma(v)
                end if
            end do
        end do
    end subroutine sssp_weighted

    function reachability(g) result(r)
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::r(:,:)
        type(geodist_result)::gd
        integer::i,j,n
        n=size(g,1)
        gd=geodist(g,.true.,.false.)
        allocate(r(n,n))
        r=0.0_dp
        do i=1,n
        do j=1,n
            if(gd%distance(i,j)<sna_inf())r(i,j)=1.0_dp
        end do
        end do
    end function reachability

    function component_dist(g,connected) result(out)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::connected
        type(component_result)::out
        character(len=16)::ct
        real(dp),allocatable::a(:,:),r(:,:)
        logical,allocatable::seen(:)
        integer::n,i,j,c,sz
        n=size(g,1)
        ct='strong'
        if(present(connected))ct=adjustl(connected)
        allocate(a(n,n))
        a=0.0_dp
        select case(trim(ct))
        case('weak')
        a=symmetrize(g,'weak')
        case('recursive')
        a=symmetrize(g,'strong')
        case('strong')
            r=reachability(g)
            do i=1,n
            do j=1,n
            if(r(i,j)/=0.0_dp.and.r(j,i)/=0.0_dp)a(i,j)=1.0_dp
            end do
            end do
        case('unilateral')
            r=reachability(g)
            do i=1,n
            do j=1,n
            if(r(i,j)/=0.0_dp.or.r(j,i)/=0.0_dp)a(i,j)=1.0_dp
            end do
            end do
        case default
        error stop 'component_dist: unknown connectedness rule'
        end select
        allocate(out%membership(n),seen(n))
        out%membership=0
        seen=.false.
        c=0
        do i=1,n
            if(seen(i))cycle
            c=c+1
            call dfs_mark(a,i,c,out%membership,seen)
        end do
        out%n_components=c
        allocate(out%csize(c))
        out%csize=0
        do i=1,n
        out%csize(out%membership(i))=out%csize(out%membership(i))+1
        end do
    end function component_dist

    recursive subroutine dfs_mark(a,v,c,memb,seen)
        real(dp),intent(in)::a(:,:)
        integer,intent(in)::v,c
        integer,intent(inout)::memb(:)
        logical,intent(inout)::seen(:)
        integer::w
        if(seen(v))return
        seen(v)=.true.
        memb(v)=c
        do w=1,size(a,1)
            if(w/=v.and.edge_present(a(v,w)).and..not.seen(w))call dfs_mark(a,w,c,memb,seen)
        end do
    end subroutine dfs_mark

    integer function components(g,connected) result(nc)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::connected
        type(component_result)::cr
        cr=component_dist(g,connected)
        nc=cr%n_components
    end function components

    function component_size_byvertex(g,connected) result(sz)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::connected
        integer,allocatable::sz(:)
        type(component_result)::cr
        integer::i
        cr=component_dist(g,connected)
        allocate(sz(size(g,1)))
        do i=1,size(sz)
        sz(i)=cr%csize(cr%membership(i))
        end do
    end function component_size_byvertex

    function component_largest_mask(g,connected) result(mask)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::connected
        logical,allocatable::mask(:)
        type(component_result)::cr
        integer::i,mx
        cr=component_dist(g,connected)
        mx=maxval(cr%csize)
        allocate(mask(size(g,1)))
        do i=1,size(mask)
        mask(i)=cr%csize(cr%membership(i))==mx
        end do
    end function component_largest_mask

    logical function is_connected(g,connected) result(ans)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::connected
        ans=components(g,connected)==1
    end function is_connected

    function isolates(g,diag) result(idx)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::diag
        integer,allocatable::idx(:)
        logical::dg
        integer::n,i,j,k
        logical,allocatable::iso(:)
        n=size(g,1)
        dg=.false.
        if(present(diag))dg=diag
        allocate(iso(n))
        iso=.true.
        do i=1,n
        do j=1,n
            if(.not.dg.and.i==j)cycle
            if(edge_present(g(i,j)).or.edge_present(g(j,i)))iso(i)=.false.
        end do
        end do
        allocate(idx(count(iso)))
        k=0
        do i=1,n
        if(iso(i))then
        k=k+1
        idx(k)=i
        end if
        end do
    end function isolates

    logical function is_isolate(g,ego,diag) result(ans)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::ego
        logical,intent(in),optional::diag
        logical::dg
        integer::j
        dg=.false.
        if(present(diag))dg=diag
        ans=.true.
        do j=1,size(g,1)
        if(.not.dg.and.j==ego)cycle
        if(edge_present(g(ego,j)).or.edge_present(g(j,ego)))then
        ans=.false.
        return
        end if
        end do
    end function is_isolate

    function cutpoints(g,mode,connected) result(cp)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::mode,connected
        logical,allocatable::cp(:)
        character(len=12)::md,ct
        real(dp),allocatable::a(:,:),sub(:,:)
        integer::n,i,j,k,base,now
        md='digraph'
        if(present(mode))md=adjustl(mode)
        ct='strong'
        if(present(connected))ct=adjustl(connected)
        n=size(g,1)
        allocate(cp(n))
        cp=.false.
        if(trim(md)=='graph')then
        a=symmetrize(g,'weak')
        base=components(a,'weak')
        else
            select case(trim(ct))
            case('weak')
            a=symmetrize(g,'weak')
            base=components(a,'weak')
            case('recursive')
            a=symmetrize(g,'strong')
            base=components(a,'weak')
            case default
            a=g
            base=components(a,'strong')
            end select
        end if
        if(n<=2)return
        do k=1,n
            allocate(sub(n-1,n-1))
            call remove_vertex(a,k,sub)
            if(trim(md)=='graph'.or.trim(ct)/='strong')then
            now=components(sub,'weak')
            else
            now=components(sub,'strong')
            end if
            cp(k)=now>base
            deallocate(sub)
        end do
    end function cutpoints

    subroutine remove_vertex(a,k,b)
        real(dp),intent(in)::a(:,:)
        integer,intent(in)::k
        real(dp),intent(out)::b(:,:)
        integer::i,j,ii,jj
        ii=0
        do i=1,size(a,1)
        if(i==k)cycle
        ii=ii+1
        jj=0
        do j=1,size(a,2)
        if(j==k)cycle
        jj=jj+1
        b(ii,jj)=a(i,j)
        end do
        end do
    end subroutine remove_vertex

    function bicomponent_dist(g,sym_rule) result(membership_count)
        ! Number of biconnected components containing each vertex.  This matrix-free
        ! interface captures the most commonly consumed part of R's bicomponent.dist.
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::sym_rule
        integer,allocatable::membership_count(:)
        real(dp),allocatable::a(:,:),sub(:,:)
        character(len=12)::sr
        integer::n,v,u,c0,c1
        sr='strong'
        if(present(sym_rule))sr=adjustl(sym_rule)
        a=symmetrize(g,sr)
        n=size(g,1)
        allocate(membership_count(n))
        membership_count=0
        ! Count blocks by edge pairs using the articulation criterion.  Each ordinary
        ! vertex receives at least one block if non-isolated; articulation vertices
        ! receive one plus the number of component increases after deletion.
        c0=components(a,'weak')
        do v=1,n
            if(.not.is_isolate(a,v))membership_count(v)=1
            if(n>2)then
                allocate(sub(n-1,n-1))
                call remove_vertex(a,v,sub)
                c1=components(sub,'weak')
                membership_count(v)=membership_count(v)+max(0,c1-c0)
                deallocate(sub)
            end if
        end do
    end function bicomponent_dist

    function kcores(g,mode,diag,cmode,ignore_eval) result(core)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::mode,cmode
        logical,intent(in),optional::diag,ignore_eval
        real(dp),allocatable::core(:),a(:,:)
        logical,allocatable::active(:)
        character(len=12)::md,cm
        logical::dg,ign,changed
        integer::n,k,i
        real(dp)::d
        n=size(g,1)
        allocate(core(n),active(n),a(n,n))
        core=0.0_dp
        active=.true.
        a=g
        md='digraph'
        if(present(mode))md=adjustl(mode)
        cm='freeman'
        if(present(cmode))cm=adjustl(cmode)
        if(trim(md)=='graph')cm='indegree'
        dg=.false.
        if(present(diag))dg=diag
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        do k=1,n
            do
                changed=.false.
                do i=1,n
                    if(.not.active(i))cycle
                    d=degree_one_active(a,i,active,cm,dg,ign)
                    if(d<real(k,dp)-sna_eps)then
                    active(i)=.false.
                    core(i)=real(k-1,dp)
                    changed=.true.
                    end if
                end do
                if(.not.changed)exit
            end do
        end do
        where(active)core=real(n,dp)
    end function kcores

    real(dp) function degree_one_active(g,i,active,cmode,diag,ignore_eval) result(d)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::i
        logical,intent(in)::active(:)
        character(*),intent(in)::cmode
        logical,intent(in)::diag,ignore_eval
        integer::j
        real(dp)::v
        d=0.0_dp
        do j=1,size(g,1)
            if(.not.active(j))cycle
            if(.not.diag.and.j==i)cycle
            select case(trim(cmode))
            case('indegree')
            if(edge_present(g(j,i)))then
            v=merge(1.0_dp,g(j,i),ignore_eval)
            d=d+v
            end if
            case('outdegree')
            if(edge_present(g(i,j)))then
            v=merge(1.0_dp,g(i,j),ignore_eval)
            d=d+v
            end if
            case default
                if(edge_present(g(j,i)))then
                v=merge(1.0_dp,g(j,i),ignore_eval)
                d=d+v
                end if
                if(edge_present(g(i,j)))then
                v=merge(1.0_dp,g(i,j),ignore_eval)
                d=d+v
                end if
            end select
        end do
    end function degree_one_active

    real(dp) function maxflow(g,src,sink,ignore_eval) result(flow)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::src,sink
        logical,intent(in),optional::ignore_eval
        real(dp),allocatable::cap(:,:),res(:,:)
        integer,allocatable::parent(:),q(:)
        logical,allocatable::seen(:)
        logical::ign,found
        integer::n,i,j,v,w,head,tail
        real(dp)::aug
        n=size(g,1)
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        allocate(cap(n,n),res(n,n),parent(n),q(n),seen(n))
        cap=0.0_dp
        do i=1,n
        do j=1,n
        if(edge_present(g(i,j)))cap(i,j)=merge(1.0_dp,max(0.0_dp,g(i,j)),ign)
        end do
        end do
        res=cap
        flow=0.0_dp
        if(src==sink)return
        do
            seen=.false.
            parent=0
            head=1
            tail=1
            q(1)=src
            seen(src)=.true.
            found=.false.
            do while(head<=tail.and..not.found)
                v=q(head)
                head=head+1
                do w=1,n
                    if(.not.seen(w).and.res(v,w)>sna_eps)then
                    seen(w)=.true.
                    parent(w)=v
                    tail=tail+1
                    q(tail)=w
                    if(w==sink)then
                    found=.true.
                    exit
                    end if
                    end if
                end do
            end do
            if(.not.found)exit
            aug=huge(1.0_dp)
            v=sink
            do while(v/=src)
            w=parent(v)
            aug=min(aug,res(w,v))
            v=w
            end do
            v=sink
            do while(v/=src)
            w=parent(v)
            res(w,v)=res(w,v)-aug
            res(v,w)=res(v,w)+aug
            v=w
            end do
            flow=flow+aug
        end do
    end function maxflow

    function neighborhood(g,order,neighborhood_type,mode,diag,thresh,partial) result(nh)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::order
        character(*),intent(in),optional::neighborhood_type,mode
        logical,intent(in),optional::diag,partial
        real(dp),intent(in),optional::thresh
        real(dp),allocatable::nh(:,:)
        character(len=12)::nt,md
        logical::dg,part
        real(dp)::thr
        integer::n,i,j
        type(geodist_result)::gd
        real(dp),allocatable::a(:,:)
        n=size(g,1)
        nt='total'
        if(present(neighborhood_type))nt=adjustl(neighborhood_type)
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        part=.true.
        if(present(partial))part=partial
        thr=0.0_dp
        if(present(thresh))thr=thresh
        allocate(a(n,n))
        a=0.0_dp
        do i=1,n
        do j=1,n
        if(.not.is_missing(g(i,j)).and.g(i,j)>thr)a(i,j)=1.0_dp
        end do
        end do
        if(trim(md)=='graph'.or.trim(nt)=='total')a=symmetrize(a,'weak')
        if(trim(nt)=='in')a=transpose(a)
        gd=geodist(a,.true.,.false.)
        allocate(nh(n,n))
        nh=0.0_dp
        do i=1,n
        do j=1,n
            if(part)then
            if(gd%distance(i,j)<=real(order,dp))nh(i,j)=1.0_dp
            else
            if(abs(gd%distance(i,j)-real(order,dp))<=sna_eps)nh(i,j)=1.0_dp
            end if
        end do
        end do
        if(.not.dg) then
            do i=1,n
                nh(i,i)=0.0_dp
            end do
        end if
    end function neighborhood

    function simmelian(g,dichotomize) result(s)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::dichotomize
        real(dp),allocatable::s(:,:)
        logical::dc
        integer::n,i,j,k
        real(dp)::v
        n=size(g,1)
        dc=.true.
        if(present(dichotomize))dc=dichotomize
        allocate(s(n,n))
        s=0.0_dp
        do i=1,n
        do j=1,n
            if(i==j)cycle
            do k=1,n
                if(k==i.or.k==j)cycle
                if(edge_present(g(i,k)).and.edge_present(g(k,j)).and.edge_present(g(j,i)) .and. &
                   edge_present(g(j,k)).and.edge_present(g(k,i)).and.edge_present(g(i,j)))then
                    if(dc)then
                    s(i,j)=1.0_dp
                    exit
                    else
                    s(i,j)=s(i,j)+1.0_dp
                    end if
                end if
            end do
        end do
        end do
    end function simmelian

    function degree(g,cmode,diag,ignore_eval,rescale) result(deg)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::diag,ignore_eval,rescale
        real(dp),allocatable::deg(:)
        character(len=12)::cm
        logical::dg,ign,rs
        integer::i
        real(dp)::ss
        cm='freeman'
        if(present(cmode))cm=adjustl(cmode)
        dg=.false.
        if(present(diag))dg=diag
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        block
            logical, allocatable :: active(:)
            allocate(active(size(g,1)))
            active=.true.
            allocate(deg(size(g,1)))
            do i=1,size(g,1)
                deg(i)=degree_one_active(g,i,active,cm,dg,ign)
            end do
        end block
        if(rs)then
        ss=sum(deg)
        if(abs(ss)>sna_eps)deg=deg/ss
        end if
    end function degree

    function dyad_census(g) result(dc)
        real(dp),intent(in)::g(:,:)
        real(dp)::dc(3)
        integer::n,i,j
        logical::a,b
        n=size(g,1)
        dc=0.0_dp
        do i=1,n-1
        do j=i+1,n
            if(is_missing(g(i,j)).or.is_missing(g(j,i)))cycle
            a=edge_present(g(i,j))
            b=edge_present(g(j,i))
            if(a.and.b)then
            dc(1)=dc(1)+1.0_dp
            else if(a.or.b)then
            dc(2)=dc(2)+1.0_dp
            else
            dc(3)=dc(3)+1.0_dp
            end if
        end do
        end do
    end function dyad_census

    real(dp) function gden(g,diag,mode,ignore_eval) result(den)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::diag,ignore_eval
        character(*),intent(in),optional::mode
        logical::dg,ign
        character(len=12)::md
        integer::n,i,j,avail
        real(dp)::cnt
        n=size(g,1)
        dg=.false.
        if(present(diag))dg=diag
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        md='digraph'
        if(present(mode))md=adjustl(mode)
        cnt=0.0_dp
        avail=0
        if(trim(md)=='graph')then
            do i=1,n
            do j=1,i
                if(.not.dg.and.i==j)cycle
                if(is_missing(g(i,j)))cycle
                avail=avail+1
                if(edge_present(g(i,j)))cnt=cnt+merge(1.0_dp,g(i,j),ign)
            end do
            end do
        else
            do i=1,n
            do j=1,n
                if(.not.dg.and.i==j)cycle
                if(is_missing(g(i,j)))cycle
                avail=avail+1
                if(edge_present(g(i,j)))cnt=cnt+merge(1.0_dp,g(i,j),ign)
            end do
            end do
        end if
        if(avail==0)then
        den=sna_nan()
        else
        den=cnt/real(avail,dp)
        end if
    end function gden

    real(dp) function connectedness(g) result(con)
        real(dp),intent(in)::g(:,:)
        type(component_result)::cr
        integer::n,i
        n=size(g,1)
        if(n<=1)then
        con=1.0_dp
        return
        end if
        cr=component_dist(g,'weak')
        con=0.0_dp
        do i=1,cr%n_components
        con=con+real(cr%csize(i)*(cr%csize(i)-1),dp)
        end do
        con=con/real(n*(n-1),dp)
    end function connectedness

    real(dp) function efficiency(g,diag) result(eff)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::diag
        logical::dg
        type(component_result)::cr
        integer::i,j,n
        real(dp)::req,maxv,edgec
        dg=.false.
        if(present(diag))dg=diag
        n=size(g,1)
        cr=component_dist(g,'weak')
        req=0.0_dp
        maxv=0.0_dp
        do i=1,cr%n_components
            req=req+real(cr%csize(i)-1,dp)
            maxv=maxv+real(cr%csize(i)*(cr%csize(i)-merge(0,1,dg))-(cr%csize(i)-1),dp)
        end do
        edgec=0.0_dp
        do i=1,n
        do j=1,n
        if(.not.dg.and.i==j)cycle
        if(edge_present(g(i,j)))edgec=edgec+g(i,j)
        end do
        end do
        if(abs(maxv)<=sna_eps)then
        eff=1.0_dp
        else
        eff=1.0_dp-(edgec-req)/maxv
        end if
    end function efficiency

    real(dp) function grecip(g,measure) result(r)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::measure
        character(len=20)::meas
        real(dp)::dc(3),den
        meas='dyadic'
        if(present(measure))meas=adjustl(measure)
        dc=dyad_census(g)
        select case(trim(meas))
        case('dyadic')
        den=sum(dc)
        r=merge((dc(1)+dc(3))/den,sna_nan(),den>0.0_dp)
        case('dyadic.nonnull')
        den=dc(1)+dc(2)
        r=merge(dc(1)/den,sna_nan(),den>0.0_dp)
        case('edgewise')
        den=2.0_dp*dc(1)+dc(2)
        r=merge(2.0_dp*dc(1)/den,sna_nan(),den>0.0_dp)
        case('edgewise.lrr')
        r=log(dc(1)*sum(dc)/(dc(1)+dc(2)/2.0_dp)**2)
        case('correlation')
        r=adjacency_recip_corr(g)
        case default
        error stop 'grecip: unknown measure'
        end select
    end function grecip

    real(dp) function adjacency_recip_corr(g) result(r)
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::x(:),y(:)
        integer::n,i,j,k
        real(dp)::mx,my,sx,sy
        n=size(g,1)
        allocate(x(n*(n-1)),y(n*(n-1)))
        k=0
        do i=1,n
        do j=1,n
        if(i==j)cycle
        if(is_missing(g(i,j)).or.is_missing(g(j,i)))cycle
        k=k+1
        x(k)=g(i,j)
        y(k)=g(j,i)
        end do
        end do
        if(k<2)then
        r=sna_nan()
        return
        end if
        mx=sum(x(:k))/k
        my=sum(y(:k))/k
        sx=sum((x(:k)-mx)**2)
        sy=sum((y(:k)-my)**2)
        if(sx<=sna_eps.or.sy<=sna_eps)then
        r=1.0_dp
        else
        r=sum((x(:k)-mx)*(y(:k)-my))/sqrt(sx*sy)
        end if
    end function adjacency_recip_corr

    real(dp) function mutuality(g) result(m)
        real(dp),intent(in)::g(:,:)
        real(dp) :: dc(3)
        dc=dyad_census(g)
        m=dc(1)
    end function mutuality

    integer function triad_classify(g,i,j,k,directed) result(tc)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::i,j,k
        logical,intent(in),optional::directed
        logical::dir
        integer::sij,sji,sjk,skj,sik,ski,m,a,n,di,dj,dk
        dir=.true.
        if(present(directed))dir=directed
        if(any([is_missing(g(i,j)),is_missing(g(j,i)),is_missing(g(j,k)),is_missing(g(k,j)),is_missing(g(i,k)), &
            & is_missing(g(k,i))]))then
        tc=-1
        return
        end if
        sij=merge(1,0,edge_present(g(i,j)))
        sji=merge(1,0,edge_present(g(j,i)))
        sjk=merge(1,0,edge_present(g(j,k)))
        skj=merge(1,0,edge_present(g(k,j)))
        sik=merge(1,0,edge_present(g(i,k)))
        ski=merge(1,0,edge_present(g(k,i)))
        if(.not.dir)then
        tc=sij+sjk+sik
        return
        end if
        m=sij*sji+sjk*skj+sik*ski
        n=(1-sij)*(1-sji)+(1-sjk)*(1-skj)+(1-sik)*(1-ski)
        a=3-m-n
        if(n==3)then
        tc=0
        else if(a==1.and.n==2)then
        tc=1
        else if(m==1.and.n==2)then
        tc=2
        else if(a==2.and.n==1)then
            di=sij+sik
            dj=sji+sjk
            dk=ski+skj
            if(di==2.or.dj==2.or.dk==2)then
            tc=3
            return
            end if
            di=sji+ski
            dj=sij+skj
            dk=sik+sjk
            if(di==2.or.dj==2.or.dk==2)then
            tc=4
            return
            end if
            tc=5
        else if(m==1.and.n==1)then
            di=sji+ski
            dj=sij+skj
            if(di==0.or.di==2.or.dj==0.or.dj==2)then
            tc=6
            else
            tc=7
            end if
        else if(a==3)then
            di=sji+ski
            dj=sij+skj
            if(di==0.or.di==2.or.dj==0.or.dj==2)then
            tc=8
            else
            tc=9
            end if
        else if(m==2.and.n==1)then
        tc=10
        else if(m==1.and.a==2)then
            di=sji+ski
            dj=sij+skj
            dk=sik+sjk
            if(di==0.or.dj==0.or.dk==0)then
            tc=11
            return
            end if
            di=sij+sik
            dj=sji+sjk
            dk=ski+skj
            if(di==0.or.dj==0.or.dk==0)then
            tc=12
            return
            end if
            tc=13
        else if(m==2.and.a==1)then
        tc=14
        else
        tc=15
        end if
    end function triad_classify

    function triad_census(g,directed) result(census)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::directed
        real(dp),allocatable::census(:)
        logical::dir
        integer::n,i,j,k,tc
        dir=.true.
        if(present(directed))dir=directed
        n=size(g,1)
        allocate(census(merge(16,4,dir)))
        census=0.0_dp
        do i=1,n-2
        do j=i+1,n-1
        do k=j+1,n
        tc=triad_classify(g,i,j,k,dir)
        if(tc>=0)census(tc+1)=census(tc+1)+1.0_dp
        end do
        end do
        end do
    end function triad_census

    real(dp) function gtrans(g,measure) result(t)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::measure
        character(len=16)::meas
        integer::n,i,j,k,risk,ok
        real(dp)::a,b,c
        meas='weak'
        if(present(measure))meas=adjustl(measure)
        n=size(g,1)
        risk=0
        ok=0
        select case(trim(meas))
        case('weak','weakcensus')
            do i=1,n
            do j=1,n
            if(j==i.or..not.edge_present(g(i,j)))cycle
            do k=1,n
            if(k==i.or.k==j.or..not.edge_present(g(j,k)))cycle
            if(is_missing(g(i,k)))cycle
            risk=risk+1
            if(edge_present(g(i,k)))ok=ok+1
            end do
            end do
            end do
        case('strong','strongcensus')
            do i=1,n
            do j=1,n
            if(j==i)cycle
            do k=1,n
            if(k==i.or.k==j)cycle
            if(is_missing(g(i,j)).or.is_missing(g(j,k)).or.is_missing(g(i,k)))cycle
            risk=risk+1
            if((edge_present(g(i,j)).and.edge_present(g(j,k)).and.edge_present(g(i,k))).or.((.not.(edge_present(g(i, &
                & j)).and.edge_present(g(j,k)))).and..not.edge_present(g(i,k))))ok=ok+1
            end do
            end do
            end do
        case('rank')
            do i=1,n
            do j=1,n
            if(j==i.or..not.edge_present(g(i,j)))cycle
            do k=1,n
            if(k==i.or.k==j.or..not.edge_present(g(j,k)).or.is_missing(g(i,k)))cycle
            risk=risk+1
            if(g(i,k)>=min(g(i,j),g(j,k)))ok=ok+1
            end do
            end do
            end do
        case('correlation')
            t=trans_corr(g)
            return
        case default
        error stop 'gtrans: unknown measure'
        end select
        if(index(trim(meas),'census')>0)then
            t=real(ok,dp)
        else if(risk==0)then
            t=1.0_dp
        else
            t=real(ok,dp)/real(risk,dp)
        end if
    end function gtrans

    real(dp) function trans_corr(g) result(r)
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::x(:),y(:)
        integer::n,i,j,k,m
        real(dp)::mx,my,sx,sy
        n=size(g,1)
        allocate(x(n*(n-1)*(max(0,n-2))))
        allocate(y(size(x)))
        m=0
        do i=1,n
        do j=1,n
        if(j==i)cycle
        do k=1,n
        if(k==i.or.k==j)cycle
        if(is_missing(g(i,k)).or.is_missing(g(i,j)).or.is_missing(g(j,k)))cycle
        m=m+1
        x(m)=g(i,k)
        y(m)=g(i,j)*g(j,k)
        end do
        end do
        end do
        if(m<2)then
        r=sna_nan()
        return
        end if
        mx=sum(x(:m))/m
        my=sum(y(:m))/m
        sx=sum((x(:m)-mx)**2)
        sy=sum((y(:m)-my)**2)
        if(sx<=sna_eps.or.sy<=sna_eps)then
        r=merge(1.0_dp,0.0_dp,all(abs(x(:m)-y(:m))<=sna_eps))
        else
        r=sum((x(:m)-mx)*(y(:m)-my))/sqrt(sx*sy)
        end if
    end function trans_corr

    real(dp) function hierarchy(g,measure) result(h)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::measure
        character(len=16)::meas
        real(dp)::rec
        meas='reciprocity'
        if(present(measure))meas=adjustl(measure)
        select case(trim(meas))
        case('reciprocity')
        rec=grecip(g,'edgewise')
        h=1.0_dp-rec
        case('krackhardt')
        h=krack_hierarchy(g)
        case default
        error stop 'hierarchy: unknown measure'
        end select
    end function hierarchy

    real(dp) function krack_hierarchy(g) result(h)
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::r(:,:)
        integer::i,j,n,asym,tot
        r=reachability(g)
        n=size(g,1)
        asym=0
        tot=0
        do i=1,n
        do j=1,n
        if(i==j)cycle
        if(r(i,j)/=0.0_dp.or.r(j,i)/=0.0_dp)then
        tot=tot+1
        if(r(i,j)/=r(j,i))asym=asym+1
        end if
        end do
        end do
        if(tot==0)then
        h=1.0_dp
        else
        h=real(asym,dp)/real(tot,dp)
        end if
    end function krack_hierarchy

    real(dp) function lubness(g) result(lub)
        ! Krackhardt upper semilattice measure: fraction of reachable unordered
        ! vertex pairs having a unique least upper bound in the reachability poset.
        real(dp),intent(in)::g(:,:)
        real(dp),allocatable::r(:,:)
        integer::n,i,j,k,l,cand,viol,tot
        logical::is_lub
        r=reachability(g)
        n=size(g,1)
        viol=0
        tot=0
        do i=1,n-1
        do j=i+1,n
        tot=tot+1
        cand=0
            do k=1,n
                if(r(i,k)==0.0_dp.or.r(j,k)==0.0_dp)cycle
                is_lub=.true.
                do l=1,n
                    if(l==k)cycle
                    if(r(i,l)/=0.0_dp.and.r(j,l)/=0.0_dp.and.r(l,k)/=0.0_dp.and.r(k,l)==0.0_dp)then
                    is_lub=.false.
                    exit
                    end if
                end do
                if(is_lub)cand=cand+1
            end do
            if(cand/=1)viol=viol+1
        end do
        end do
        if(tot==0)then
        lub=1.0_dp
        else
        lub=1.0_dp-real(viol,dp)/real(tot,dp)
        end if
    end function lubness

    function structure_statistics(g) result(stat)
        real(dp),intent(in)::g(:,:)
        real(dp)::stat(8)
        stat(1)=gden(g,.false.,'digraph',.true.)
        stat(2)=grecip(g,'edgewise')
        stat(3)=gtrans(g,'weak')
        stat(4)=connectedness(g)
        stat(5)=hierarchy(g,'krackhardt')
        stat(6)=lubness(g)
        stat(7)=efficiency(g,.false.)
        stat(8)=real(components(g,'strong'),dp)
    end function structure_statistics

    function kpath_census(g,maxlen,directed,by_vertex,copaths,dyadpaths) result(out)
        real(dp),intent(in)::g(:,:)
        integer,intent(in),optional::maxlen
        logical,intent(in),optional::directed,by_vertex
        integer,intent(in),optional::copaths,dyadpaths
        type(path_census_result)::out
        integer::n,ml,src
        logical::dir,bv
        integer::cp,dpd
        logical,allocatable::used(:)
        integer,allocatable::path(:)
        n=size(g,1)
        ml=min(3,n-1)
        if(present(maxlen))ml=min(maxlen,n-1)
        dir=.true.
        if(present(directed))dir=directed
        bv=.true.
        if(present(by_vertex))bv=by_vertex
        cp=0
        if(present(copaths))cp=copaths
        dpd=0
        if(present(dyadpaths))dpd=dyadpaths
        allocate(out%count(ml))
        out%count=0.0_dp
        if(bv)then
        allocate(out%vertex_count(ml,n))
        out%vertex_count=0.0_dp
        else
        allocate(out%vertex_count(0,0))
        end if
        if(dpd>0)then
        allocate(out%dyad_count(ml,n,n))
        out%dyad_count=0.0_dp
        else
        allocate(out%dyad_count(0,0,0))
        end if
        allocate(used(n),path(n))
        do src=1,n
        used=.false.
        used(src)=.true.
        path(1)=src
        call path_dfs(g,src,src,0,ml,dir,used,path,out,bv,dpd)
        end do
    end function kpath_census

    recursive subroutine path_dfs(g,src,v,len,maxlen,dir,used,path,out,bv,dpd)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::src,v,len,maxlen
        logical,intent(in)::dir,bv
        logical,intent(inout)::used(:)
        integer,intent(inout)::path(:)
        type(path_census_result),intent(inout)::out
        integer,intent(in)::dpd
        integer::w,l,u
        if(len>=maxlen)return
        do w=1,size(g,1)
            if(w==v.or.used(w))cycle
            if(dir)then
            if(.not.edge_present(g(v,w)))cycle
            else
            if(.not.(edge_present(g(v,w)).or.edge_present(g(w,v))))cycle
            end if
            l=len+1
            used(w)=.true.
            path(l+1)=w
            out%count(l)=out%count(l)+1.0_dp
            if(bv) then
                do u=1,l+1
                    out%vertex_count(l,path(u))=out%vertex_count(l,path(u))+1.0_dp
                end do
            end if
            if(dpd>0)out%dyad_count(l,src,w)=out%dyad_count(l,src,w)+1.0_dp
            call path_dfs(g,src,w,l,maxlen,dir,used,path,out,bv,dpd)
            used(w)=.false.
        end do
    end subroutine path_dfs

    function kcycle_census(g,maxlen,directed) result(count)
        real(dp),intent(in)::g(:,:)
        integer,intent(in),optional::maxlen
        logical,intent(in),optional::directed
        real(dp),allocatable::count(:)
        integer::n,ml,s
        logical::dir
        logical,allocatable::used(:)
        n=size(g,1)
        ml=min(3,n)
        if(present(maxlen))ml=min(maxlen,n)
        dir=.true.
        if(present(directed))dir=directed
        allocate(count(max(0,ml-1)))
        count=0.0_dp
        allocate(used(n))
        do s=1,n
        used=.false.
        used(s)=.true.
        call cycle_dfs(g,s,s,0,ml,dir,used,count)
        end do
        ! Each directed cycle is encountered once per starting vertex; each undirected
        ! cycle also in two orientations.
        do s=2,ml
            if(dir)then
            count(s-1)=count(s-1)/real(s,dp)
            else
            count(s-1)=count(s-1)/real(2*s,dp)
            end if
        end do
    end function kcycle_census

    recursive subroutine cycle_dfs(g,src,v,len,maxlen,dir,used,count)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::src,v,len,maxlen
        logical,intent(in)::dir
        logical,intent(inout)::used(:)
        real(dp),intent(inout)::count(:)
        integer::w,l
        if(len>=maxlen)return
        do w=1,size(g,1)
            if(dir)then
            if(.not.edge_present(g(v,w)))cycle
            else
            if(.not.(edge_present(g(v,w)).or.edge_present(g(w,v))))cycle
            end if
            l=len+1
            if(w==src.and.l>=2)then
            if(l<=maxlen)count(l-1)=count(l-1)+1.0_dp
            cycle
            end if
            if(used(w))cycle
            used(w)=.true.
            call cycle_dfs(g,src,w,l,maxlen,dir,used,count)
            used(w)=.false.
        end do
    end subroutine cycle_dfs

    function clique_census(g) result(count_by_size)
        real(dp),intent(in)::g(:,:)
        integer,allocatable::count_by_size(:)
        integer::n
        logical,allocatable::chosen(:)
        n=size(g,1)
        allocate(count_by_size(n))
        count_by_size=0
        allocate(chosen(n))
        chosen=.false.
        call clique_enum(g,1,0,chosen,count_by_size)
    end function clique_census

    recursive subroutine clique_enum(g,start,sz,chosen,count)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::start,sz
        logical,intent(inout)::chosen(:)
        integer,intent(inout)::count(:)
        integer::v,u
        logical::ok
        do v=start,size(g,1)
            ok=.true.
            do u=1,size(g,1)
            if(chosen(u).and..not.(edge_present(g(v,u)).or.edge_present(g(u,v))))then
            ok=.false.
            exit
            end if
            end do
            if(.not.ok)cycle
            chosen(v)=.true.
            count(sz+1)=count(sz+1)+1
            call clique_enum(g,v+1,sz+1,chosen,count)
            chosen(v)=.false.
        end do
    end subroutine clique_enum

end module sna_graph
