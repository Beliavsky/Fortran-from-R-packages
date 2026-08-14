! Node-level indices translated from R/sna nli.R and src/nli.c.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_centrality
    use sna_kinds, only : dp, sna_eps, sna_inf, is_missing
    use sna_types, only : geodist_result
    use sna_graph, only : edge_present, geodist, maxflow, reachability, is_isolate, degree
    use sna_prep, only : symmetrize, make_stochastic
    use sna_linalg, only : inverse_matrix, dominant_eigenvector
    implicit none
    private

    public :: betweenness, closeness, bonpow, evcent, graphcent, gilschmidt
    public :: stresscent, loadcent, flowbet, infocent, prestige
    public :: centralization_from_scores

contains

    function betweenness(g,cmode,ignore_eval,rescale) result(bet)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::bet(:),a(:,:)
        character(len=16)::cm
        logical::ign,rs,undir
        real(dp)::s
        cm='directed'
        if(present(cmode))cm=adjustl(cmode)
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        undir=trim(cm)=='undirected'
        if(undir)then
        a=symmetrize(g,'weak')
        else
        a=g
        end if
        select case(trim(cm))
        case('directed','undirected')
        call brandes_standard(a,ign,bet)
        case('endpoints')
        call brandes_endpoints(a,ign,bet)
        case('proximalsrc')
        call bet_proximal_source(a,ign,bet)
        case('proximaltar')
        call bet_proximal_target(a,ign,bet)
        case('proximalsum')
            call bet_proximal_source(a,ign,bet)
            block
            real(dp),allocatable::z(:)
            call bet_proximal_target(a,ign,z)
            bet=bet+z
            end block
        case('lengthscaled')
        call bet_lengthscaled(a,ign,bet)
        case('linearscaled')
        call bet_linearscaled(a,ign,bet)
        case default
        error stop 'betweenness: unknown cmode'
        end select
        if(undir)bet=bet/2.0_dp
        if(rs)then
        s=sum(bet)
        if(abs(s)>sna_eps)bet=bet/s
        end if
    end function betweenness

    subroutine brandes_standard(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        integer::n,s,v,w,top
        real(dp),allocatable::dist(:),sigma(:),delta(:)
        integer,allocatable::stack(:),pred(:,:),npred(:)
        n=size(g,1)
        allocate(bet(n),dist(n),sigma(n),delta(n),stack(n),pred(n,n),npred(n))
        bet=0.0_dp
        do s=1,n
            call shortest_with_pred(g,s,ignore_eval,dist,sigma,stack,top,pred,npred)
            delta=0.0_dp
            do while(top>0)
                w=stack(top)
                top=top-1
                do v=1,npred(w)
                if(sigma(w)>0.0_dp)delta(pred(v,w))=delta(pred(v,w))+sigma(pred(v,w))/sigma(w)*(1.0_dp+delta(w))
                end do
                if(w/=s)bet(w)=bet(w)+delta(w)
            end do
        end do
    end subroutine brandes_standard

    subroutine brandes_endpoints(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        integer::n,s,v,w,top,reachable
        real(dp),allocatable::dist(:),sigma(:),delta(:)
        integer,allocatable::stack(:),pred(:,:),npred(:)
        n=size(g,1)
        allocate(bet(n),dist(n),sigma(n),delta(n),stack(n),pred(n,n),npred(n))
        bet=0.0_dp
        do s=1,n
            call shortest_with_pred(g,s,ignore_eval,dist,sigma,stack,top,pred,npred)
            reachable=top-1
            bet(s)=bet(s)+real(max(0,reachable),dp)
            delta=0.0_dp
            do while(top>0)
            w=stack(top)
            top=top-1
            do v=1,npred(w)
            if(sigma(w)>0.0_dp)delta(pred(v,w))=delta(pred(v,w))+sigma(pred(v,w))/sigma(w)*(1.0_dp+delta(w))
            end do
            if(w/=s)bet(w)=bet(w)+delta(w)+1.0_dp
            end do
        end do
    end subroutine brandes_endpoints

    subroutine shortest_with_pred(g,s,ignore_eval,dist,sigma,stack,top,pred,npred)
        real(dp),intent(in)::g(:,:)
        integer,intent(in)::s
        logical,intent(in)::ignore_eval
        real(dp),intent(out)::dist(:),sigma(:)
        integer,intent(out)::stack(:),top,pred(:,:),npred(:)
        integer::n,v,w,i,bestv,qh,qt
        integer,allocatable::q(:)
        logical,allocatable::done(:)
        real(dp)::alt,best,wt
        n=size(g,1)
        dist=sna_inf()
        sigma=0.0_dp
        npred=0
        pred=0
        top=0
        dist(s)=0.0_dp
        sigma(s)=1.0_dp
        if(ignore_eval)then
            allocate(q(n))
            qh=1
            qt=1
            q(1)=s
            do while(qh<=qt)
            v=q(qh)
            qh=qh+1
            top=top+1
            stack(top)=v
                do w=1,n
                if(w==v.or..not.edge_present(g(v,w)))cycle
                    if(dist(w)==sna_inf())then
                    dist(w)=dist(v)+1.0_dp
                    qt=qt+1
                    q(qt)=w
                    end if
                    if(abs(dist(w)-(dist(v)+1.0_dp))<=sna_eps)then
                    sigma(w)=sigma(w)+sigma(v)
                    npred(w)=npred(w)+1
                    pred(npred(w),w)=v
                    end if
                end do
            end do
        else
            allocate(done(n))
            done=.false.
            do i=1,n
                best=sna_inf()
                bestv=0
                do v=1,n
                if(.not.done(v).and.dist(v)<best)then
                best=dist(v)
                bestv=v
                end if
                end do
                if(bestv==0)exit
                v=bestv
                done(v)=.true.
                top=top+1
                stack(top)=v
                do w=1,n
                if(w==v.or..not.edge_present(g(v,w)))cycle
                wt=g(v,w)
                if(wt<0.0_dp)error stop 'negative edge value in centrality'
                alt=dist(v)+wt
                    if(alt<dist(w)-sna_eps)then
                    dist(w)=alt
                    sigma(w)=sigma(v)
                    npred(w)=1
                    pred(1,w)=v
                    else if(abs(alt-dist(w))<=sna_eps)then
                    sigma(w)=sigma(w)+sigma(v)
                    npred(w)=npred(w)+1
                    pred(npred(w),w)=v
                    end if
                end do
            end do
        end if
    end subroutine shortest_with_pred

    subroutine bet_proximal_source(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        integer::n,s,v,w,top,p
        real(dp),allocatable::d(:),sig(:)
        integer,allocatable::st(:),pr(:,:),np(:)
        n=size(g,1)
        allocate(bet(n),d(n),sig(n),st(n),pr(n,n),np(n))
        bet=0.0_dp
        do s=1,n
        call shortest_with_pred(g,s,ignore_eval,d,sig,st,top,pr,np)
        do w=1,n
        do p=1,np(w)
        v=pr(p,w)
        if(v/=s.and.sig(w)>0)bet(v)=bet(v)+sig(v)/sig(w)
        end do
        end do
        end do
    end subroutine bet_proximal_source

    subroutine bet_proximal_target(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        ! Reverse-graph proximal source is equivalent to proximal target.
        call bet_proximal_source(transpose(g),ignore_eval,bet)
    end subroutine bet_proximal_target

    subroutine bet_lengthscaled(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        integer::n,s,v,w,top,p
        real(dp),allocatable::d(:),sig(:),del(:)
        integer,allocatable::st(:),pr(:,:),np(:)
        n=size(g,1)
        allocate(bet(n),d(n),sig(n),del(n),st(n),pr(n,n),np(n))
        bet=0
        do s=1,n
        call shortest_with_pred(g,s,ignore_eval,d,sig,st,top,pr,np)
        del=0
            do while(top>0)
            w=st(top)
            top=top-1
            if(d(w)>0.and.d(w)<sna_inf())then
            do p=1,np(w)
            v=pr(p,w)
            if(sig(w)>0)del(v)=del(v)+sig(v)/sig(w)*(1.0_dp/d(w)+del(w))
            end do
            end if
            if(w/=s)bet(w)=bet(w)+del(w)
            end do
        end do
    end subroutine bet_lengthscaled

    subroutine bet_linearscaled(g,ignore_eval,bet)
        real(dp),intent(in)::g(:,:)
        logical,intent(in)::ignore_eval
        real(dp),allocatable,intent(out)::bet(:)
        ! Linear scaling in sna weights each source-target contribution by reciprocal
        ! geodesic length times standard dependency; compute explicitly over triples.
        integer::n,s,t,v
        type(geodist_result)::gd
        n=size(g,1)
        allocate(bet(n))
        bet=0.0_dp
        gd=geodist(g,ignore_eval,.true.)
        do s=1,n
        do t=1,n
        if(s==t.or.gd%distance(s,t)>=sna_inf().or.gd%distance(s,t)<=0)cycle
        do v=1,n
        if(v==s.or.v==t)cycle
            if(abs(gd%distance(s,v)+gd%distance(v,t)-gd%distance(s,t))<=sna_eps.and.gd%counts(s,t)>0)bet(v)=bet(v)+ &
                & (gd%counts(s,v)*gd%counts(v,t)/gd%counts(s,t))/gd%distance(s,t)
        end do
        end do
        end do
    end subroutine bet_linearscaled

    function closeness(g,cmode,ignore_eval,rescale) result(c)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::c(:),a(:,:)
        character(len=16)::cm
        logical::ign,rs
        type(geodist_result)::gd
        integer::n,i,j,r
        real(dp)::s,ids
        cm='directed'
        if(present(cmode))cm=adjustl(cmode)
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        if(trim(cm)=='undirected'.or.trim(cm)=='suminvundir')then
        a=symmetrize(g,'weak')
        else
        a=g
        end if
        n=size(a,1)
        gd=geodist(a,ign,.false.)
        allocate(c(n))
        c=0
        do i=1,n
            select case(trim(cm))
            case('directed','undirected')
            s=0
            do j=1,n
            if(j/=i.and.gd%distance(i,j)<sna_inf())s=s+gd%distance(i,j)
            end do
            if(s>0)c(i)=real(n-1,dp)/s
            case('suminvdir','suminvundir')
            ids=0
            do j=1,n
            if(j/=i.and.gd%distance(i,j)>0.and.gd%distance(i,j)<sna_inf())ids=ids+1.0_dp/gd%distance(i,j)
            end do
            c(i)=ids/real(max(1,n-1),dp)
            case('gil-schmidt')
            ids=0
            r=0
            do j=1,n
            if(j/=i.and.gd%distance(i,j)>0.and.gd%distance(i,j)<sna_inf())then
            ids=ids+1.0_dp/gd%distance(i,j)
            r=r+1
            end if
            end do
            if(r>0)c(i)=ids/real(r,dp)
            case default
            error stop 'closeness: unknown cmode'
            end select
        end do
        if(rs)then
        s=sum(c)
        if(abs(s)>sna_eps)c=c/s
        end if
    end function closeness

    function bonpow(g,exponent,diag,rescale,tol) result(ev)
        real(dp),intent(in)::g(:,:)
        real(dp),intent(in),optional::exponent,tol
        logical,intent(in),optional::diag,rescale
        real(dp),allocatable::ev(:),a(:,:),m(:,:),inv(:,:),ones(:)
        real(dp)::ex,nrm,s
        logical::dg,rs
        integer::n,i,info
        n=size(g,1)
        ex=1.0_dp
        if(present(exponent))ex=exponent
        dg=.false.
        if(present(diag))dg=diag
        rs=.false.
        if(present(rescale))rs=rescale
        allocate(a(n,n),m(n,n),inv(n,n),ev(n),ones(n))
        a=g
        if(.not.dg) then
            do i=1,n
                a(i,i)=0.0_dp
            end do
        end if
        m=-ex*a
        do i=1,n
            m(i,i)=m(i,i)+1.0_dp
        end do
        block
            real(dp) :: at
            at=1.0e-7_dp
            if(present(tol)) at=tol
            call inverse_matrix(m,inv,info,at)
        end block
        if(info/=0)then
        ev=0
        return
        end if
        ones=1.0_dp
        ev=matmul(matmul(inv,a),ones)
        nrm=sqrt(sum(ev*ev))
        if(nrm>0)ev=ev*sqrt(real(n,dp))/nrm
        if(rs)then
        s=sum(ev)
        if(abs(s)>sna_eps)ev=ev/s
        end if
    end function bonpow

    function evcent(g,diag,rescale,ignore_eval,tol,maxiter) result(ev)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::diag,rescale,ignore_eval
        real(dp),intent(in),optional::tol
        integer,intent(in),optional::maxiter
        real(dp),allocatable::ev(:),a(:,:)
        real(dp)::val,s
        logical::dg,rs,ign,ok
        integer::n,i
        n=size(g,1)
        allocate(a(n,n),ev(n))
        a=g
        dg=.false.
        if(present(diag))dg=diag
        rs=.false.
        if(present(rescale))rs=rescale
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        if(.not.dg) then
            do i=1,n
                a(i,i)=0.0_dp
            end do
        end if
        if(ign) then
            where(.not.is_missing(a).and.a/=0.0_dp) a=1.0_dp
        end if
        where(is_missing(a)) a=0.0_dp
        call dominant_eigenvector(a,ev,val,tol,maxiter,ok)
        if(rs)then
        s=sum(ev)
        if(abs(s)>sna_eps)ev=ev/s
        end if
    end function evcent

    function graphcent(g,cmode,ignore_eval,rescale) result(gc)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::gc(:),a(:,:)
        character(len=16)::cm
        logical::ign,rs
        type(geodist_result)::gd
        integer::n,i
        real(dp)::mx,s
        cm='directed'
        if(present(cmode))cm=adjustl(cmode)
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        if(trim(cm)=='undirected')then
        a=symmetrize(g,'weak')
        else
        a=g
        end if
        n=size(a,1)
        gd=geodist(a,ign,.false.)
        allocate(gc(n))
        do i=1,n
        mx=maxval(gd%distance(i,:))
        if(mx>=sna_inf())then
        gc(i)=0
        else if(mx>0)then
        gc(i)=1.0_dp/mx
        else
        gc(i)=0
        end if
        end do
        if(rs)then
        s=sum(gc)
        if(abs(s)>sna_eps)gc=gc/s
        end if
    end function graphcent

    function gilschmidt(g,normalize) result(gs)
        real(dp),intent(in)::g(:,:)
        logical,intent(in),optional::normalize
        real(dp),allocatable::gs(:)
        logical::norm
        type(geodist_result)::gd
        integer::n,i,j,r
        norm=.true.
        if(present(normalize))norm=normalize
        n=size(g,1)
        gd=geodist(g,.true.,.false.)
        allocate(gs(n))
        gs=0
        do i=1,n
        r=0
        do j=1,n
        if(j/=i.and.gd%distance(i,j)>0.and.gd%distance(i,j)<sna_inf())then
        gs(i)=gs(i)+1.0_dp/gd%distance(i,j)
        r=r+1
        end if
        end do
        if(norm.and.r>0)gs(i)=gs(i)/real(r,dp)
        end do
    end function gilschmidt

    function stresscent(g,cmode,ignore_eval,rescale) result(stress)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::stress(:),a(:,:)
        character(len=16)::cm
        logical::ign,rs,undir
        type(geodist_result)::gd
        integer::n,i,j,k
        real(dp)::s
        cm='directed'
        if(present(cmode))cm=adjustl(cmode)
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        undir=trim(cm)=='undirected'
        if(undir)then
        a=symmetrize(g,'weak')
        else
        a=g
        end if
        n=size(a,1)
        gd=geodist(a,ign,.true.)
        allocate(stress(n))
        stress=0
        do i=1,n
        do j=1,n
        do k=1,n
            if(i==j.or.i==k.or.j==k)cycle
            if(gd%distance(j,k)<sna_inf().and.abs(gd%distance(j,i)+gd%distance(i,k)-gd%distance(j,k))<=sna_eps)stress(i)= &
                & stress(i)+gd%counts(j,i)*gd%counts(i,k)
        end do
        end do
        end do
        if(undir)stress=stress/2
        if(rs)then
        s=sum(stress)
        if(abs(s)>sna_eps)stress=stress/s
        end if
    end function stresscent

    function loadcent(g,cmode,ignore_eval,rescale) result(load)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::load(:),a(:,:),d(:),sig(:),delta(:)
        integer,allocatable::st(:),pr(:,:),np(:)
        character(len=16)::cm
        logical::ign,rs,undir
        integer::n,s,w,v,p,top
        real(dp)::sm
        cm='directed'
        if(present(cmode))cm=adjustl(cmode)
        ign=.true.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        undir=trim(cm)=='undirected'
        if(undir)then
        a=symmetrize(g,'weak')
        else
        a=transpose(g)
        end if
        n=size(a,1)
        allocate(load(n),d(n),sig(n),delta(n),st(n),pr(n,n),np(n))
        load=0
        do s=1,n
        call shortest_with_pred(a,s,ign,d,sig,st,top,pr,np)
        delta=1.0_dp
            do while(top>0)
            w=st(top)
            top=top-1
            do p=1,np(w)
            v=pr(p,w)
            if(np(w)>0)delta(v)=delta(v)+delta(w)/real(np(w),dp)
            end do
            load(w)=load(w)+delta(w)
            end do
        end do
        if(undir)load=load/2
        if(rs)then
        sm=sum(load)
        if(abs(sm)>sna_eps)load=load/sm
        end if
    end function loadcent

    function flowbet(g,gmode,cmode,ignore_eval,rescale) result(flo)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::gmode,cmode
        logical,intent(in),optional::ignore_eval,rescale
        real(dp),allocatable::flo(:),mf(:,:),sub(:,:)
        character(len=12)::gm,cm
        logical::ign,rs
        integer::n,i,j,k,j2,k2
        real(dp)::red,s
        gm='digraph'
        if(present(gmode))gm=adjustl(gmode)
        cm='rawflow'
        if(present(cmode))cm=adjustl(cmode)
        ign=.false.
        if(present(ignore_eval))ign=ignore_eval
        rs=.false.
        if(present(rescale))rs=rescale
        n=size(g,1)
        allocate(mf(n,n),flo(n))
        mf=0
        flo=0
        do j=1,n
        do k=1,n
        if(j/=k)mf(j,k)=maxflow(g,j,k,ign)
        end do
        end do
        do i=1,n
        allocate(sub(n-1,n-1))
        call remove_v(g,i,sub)
        do j=1,n
        do k=1,n
            if(i==j.or.i==k.or.j==k)cycle
            if(trim(gm)=='graph'.and.j>k)cycle
            if(mf(j,k)<=sna_eps)cycle
            j2=j-merge(1,0,j>i)
            k2=k-merge(1,0,k>i)
            red=maxflow(sub,j2,k2,ign)
            select case(trim(cm))
            case('fracflow')
            flo(i)=flo(i)+(mf(j,k)-red)/mf(j,k)
            case default
            flo(i)=flo(i)+mf(j,k)-red
            end select
        end do
        end do
        deallocate(sub)
        end do
        if(trim(cm)=='normflow')then
            do i=1,n
            s=sum(mf)-sum(mf(i,:))-sum(mf(:,i))
            if(s>0)flo(i)=flo(i)/s*merge(2.0_dp,1.0_dp,trim(gm)=='graph')
            end do
        end if
        if(rs)then
        s=sum(flo)
        if(abs(s)>sna_eps)flo=flo/s
        end if
    end function flowbet

    subroutine remove_v(a,k,b)
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
    end subroutine remove_v

    function infocent(g,cmode,diag,rescale,tol) result(ic)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::diag,rescale
        real(dp),intent(in),optional::tol
        real(dp),allocatable::ic(:),m(:,:),a(:,:),ainv(:,:),sub(:,:),sinv(:,:)
        logical,allocatable::iso(:)
        integer,allocatable::ix(:)
        character(len=12)::cm
        logical::dg,rs
        integer::n,i,j,k,info,nn
        real(dp)::tr,r,s
        n=size(g,1)
        cm='weak'
        if(present(cmode))cm=adjustl(cmode)
        dg=.false.
        if(present(diag))dg=diag
        rs=.false.
        if(present(rescale))rs=rescale
        allocate(m(n,n))
        m=g
        if(maxval(abs(m-transpose(m)),mask=.not.is_missing(m))>sna_eps)m=symmetrize(m,cm)
        allocate(iso(n))
        do i=1,n
        iso(i)=is_isolate(m,i,dg)
        end do
        allocate(ix(count(.not.iso)))
        k=0
        do i=1,n
        if(.not.iso(i))then
        k=k+1
        ix(k)=i
        end if
        end do
        allocate(ic(n))
        ic=0
        if(k==0)return
        nn=k
        allocate(sub(nn,nn),a(nn,nn),sinv(nn,nn))
        do i=1,nn
        do j=1,nn
        sub(i,j)=m(ix(i),ix(j))
        end do
        end do
        a=1.0_dp-sub
        where(sub==0.0_dp)a=1.0_dp
        do i=1,nn
        a(i,i)=1.0_dp+sum(sub(i,:),mask=.not.is_missing(sub(i,:)))
        end do
        where(is_missing(a))a=1.0_dp
        block
            real(dp) :: at
            at=1.0e-20_dp
            if(present(tol)) at=tol
            call inverse_matrix(a,sinv,info,at)
        end block
        if(info/=0)return
        tr=0
        do i=1,nn
            tr=tr+sinv(i,i)
        end do
        do i=1,nn
        r=sum(sinv(i,:))
        ic(ix(i))=1.0_dp/(sinv(i,i)+(tr-2.0_dp*r)/real(nn,dp))
        end do
        if(rs)then
        s=sum(ic)
        if(abs(s)>sna_eps)ic=ic/s
        end if
    end function infocent

    function prestige(g,cmode,diag,rescale) result(p)
        real(dp),intent(in)::g(:,:)
        character(*),intent(in),optional::cmode
        logical,intent(in),optional::diag,rescale
        real(dp),allocatable::p(:),d(:,:),r(:,:)
        character(len=28)::cm
        logical::dg,rs
        integer::n,i,j
        real(dp)::s,den,num
        type(geodist_result)::gd
        cm='indegree'
        if(present(cmode))cm=adjustl(cmode)
        dg=.false.
        if(present(diag))dg=diag
        rs=.false.
        if(present(rescale))rs=rescale
        n=size(g,1)
        select case(trim(cm))
        case('indegree')
        p=degree(g,'indegree',dg,.false.,.false.)
        case('indegree.rownorm')
        d=make_stochastic(g,'row')
        p=degree(d,'indegree',dg,.false.,.false.)
        case('indegree.colnorm')
        d=make_stochastic(g,'col')
        p=degree(d,'indegree',dg,.false.,.false.)
        case('indegree.rowcolnorm')
        d=make_stochastic(g,'rowcol')
        p=degree(d,'indegree',dg,.false.,.false.)
        case('eigenvector')
        p=evcent(transpose(g),dg,.false.,.false.)
        case('eigenvector.rownorm')
        d=make_stochastic(g,'row')
        p=evcent(transpose(d),dg,.false.,.false.)
        case('eigenvector.colnorm')
        d=make_stochastic(g,'col')
        p=evcent(transpose(d),dg,.false.,.false.)
        case('eigenvector.rowcolnorm')
        d=make_stochastic(g,'rowcol')
        p=evcent(transpose(d),dg,.false.,.false.)
        case('domain')
        r=reachability(g)
        allocate(p(n))
        do j=1,n
        p(j)=sum(r(:,j))-1.0_dp
        end do
        case('domain.proximity')
        gd=geodist(g,.true.,.true.)
        allocate(p(n))
        p=0
        do j=1,n
        num=0
        den=0
        do i=1,n
        if(i/=j.and.gd%counts(i,j)>0)then
        num=num+1
        den=den+gd%distance(i,j)
        end if
        end do
        if(den>0)p(j)=num*num/(den*real(max(1,n-1),dp))
        end do
        case default
        error stop 'prestige: unknown cmode'
        end select
        if(rs)then
        s=sum(p)
        if(abs(s)>sna_eps)p=p/s
        end if
    end function prestige

    pure real(dp) function centralization_from_scores(scores,theoretical_max_deviation) result(c)
        real(dp),intent(in)::scores(:)
        real(dp),intent(in),optional::theoretical_max_deviation
        real(dp)::mx
        mx=maxval(scores)
        c=sum(mx-scores)
        if(present(theoretical_max_deviation))then
        if(abs(theoretical_max_deviation)>sna_eps)c=c/theoretical_max_deviation
        end if
    end function centralization_from_scores

end module sna_centrality
