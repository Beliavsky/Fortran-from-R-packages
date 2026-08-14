! Random graph generation translated from R/sna randomgraph.R and src/randomgraph.c.
! Upstream copyright (C) Carter T. Butts.
! Licensed under GPL-2.0-or-later; see COPYING.
module sna_random
    use sna_kinds, only : dp, sna_eps
    use sna_prep, only : symmetrize
    implicit none
    private

    public :: rgraph, rgnm, rguman, rgnmix_probability, rgnmix_exact
    public :: rewire_ud, rewire_ws, rgws, rgbn_mcmc

contains

    function rgraph(n,tprob,mode,diag) result(g)
        integer,intent(in)::n
        real(dp),intent(in),optional::tprob(:,:)
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag
        real(dp),allocatable::g(:,:)
        character(len=12)::md
        logical::dg
        real(dp)::u,p
        integer::i,j
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        allocate(g(n,n))
        g=0.0_dp
        if(trim(md)=='graph')then
            do i=1,n
            do j=i,n
            if(.not.dg.and.i==j)cycle
            p=0.5_dp
            if(present(tprob))p=tprob(i,j)
            call random_number(u)
            if(u<p)then
            g(i,j)=1
            g(j,i)=1
            end if
            end do
            end do
        else
            do i=1,n
            do j=1,n
            if(.not.dg.and.i==j)cycle
            p=0.5_dp
            if(present(tprob))p=tprob(i,j)
            call random_number(u)
            if(u<p)g(i,j)=1
            end do
            end do
        end if
    end function rgraph

    function rgnm(nv,m,mode,diag) result(g)
        integer,intent(in)::nv,m
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag
        real(dp),allocatable::g(:,:)
        character(len=12)::md
        logical::dg
        integer::maxe,k,i,j,cnt,pick
        integer,allocatable::ri(:),rj(:),ord(:)
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        allocate(g(nv,nv))
        g=0
        if(trim(md)=='graph')then
        maxe=nv*(nv-1)/2+merge(nv,0,dg)
        else
        maxe=nv*(nv-1)+merge(nv,0,dg)
        end if
        if(m<0.or.m>maxe)error stop 'rgnm: too many edges requested'
        allocate(ri(maxe),rj(maxe),ord(maxe))
        cnt=0
        if(trim(md)=='graph')then
            do i=1,nv
            do j=i,nv
            if(.not.dg.and.i==j)cycle
            cnt=cnt+1
            ri(cnt)=i
            rj(cnt)=j
            end do
            end do
        else
            do i=1,nv
            do j=1,nv
            if(.not.dg.and.i==j)cycle
            cnt=cnt+1
            ri(cnt)=i
            rj(cnt)=j
            end do
            end do
        end if
        do k=1,maxe
        ord(k)=k
        end do
        call shuffle_int(ord)
        do k=1,m
        pick=ord(k)
        i=ri(pick)
        j=rj(pick)
        g(i,j)=1
        if(trim(md)=='graph')g(j,i)=1
        end do
    end function rgnm

    function rguman(nv,mut,asym,null,method) result(g)
        integer,intent(in)::nv
        real(dp),intent(in)::mut,asym,null
        character(*),intent(in),optional::method
        real(dp),allocatable::g(:,:)
        character(len=16)::meth
        integer::nd,mc,ac,nc,i,j,k
        integer,allocatable::ri(:),rj(:),state(:),ord(:)
        real(dp)::sm,u,pm,pa
        meth='probability'
        if(present(method))meth=adjustl(method)
        nd=nv*(nv-1)/2
        allocate(g(nv,nv),ri(nd),rj(nd),state(nd),ord(nd))
        g=0
        k=0
        do i=1,nv-1
        do j=i+1,nv
        k=k+1
        ri(k)=i
        rj(k)=j
        end do
        end do
        if(trim(meth)=='exact')then
            mc=nint(mut)
            ac=nint(asym)
            nc=nint(null)
            if(mc+ac+nc/=nd)error stop 'rguman exact: dyad counts must sum to choose(nv,2)'
        else
            sm=mut+asym+null
            if(sm<=0)error stop 'rguman: probabilities sum to zero'
            pm=mut/sm
            pa=asym/sm
            mc=0
            ac=0
            nc=0
            do k=1,nd
            call random_number(u)
            if(u<pm)then
            state(k)=1
            mc=mc+1
            else if(u<pm+pa)then
            state(k)=2
            ac=ac+1
            else
            state(k)=3
            nc=nc+1
            end if
            end do
        end if
        if(trim(meth)=='exact')then
        do k=1,nd
        ord(k)=k
        end do
        call shuffle_int(ord)
        state=3
        do k=1,mc
        state(ord(k))=1
        end do
        do k=mc+1,mc+ac
        state(ord(k))=2
        end do
        end if
        do k=1,nd
        i=ri(k)
        j=rj(k)
        select case(state(k))
        case(1)
        g(i,j)=1
        g(j,i)=1
        case(2)
        call random_number(u)
        if(u<0.5_dp)then
        g(i,j)=1
        else
        g(j,i)=1
        end if
        end select
        end do
    end function rguman

    function rgnmix_probability(types,mix,mode,diag) result(g)
        integer,intent(in)::types(:)
        real(dp),intent(in)::mix(:,:)
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag
        real(dp),allocatable::g(:,:),p(:,:)
        integer::n,i,j
        n=size(types)
        allocate(p(n,n))
        do i=1,n
        do j=1,n
        p(i,j)=mix(types(i),types(j))
        end do
        end do
        g=rgraph(n,p,mode,diag)
    end function rgnmix_probability

    function rgnmix_exact(types,mix,mode,diag) result(g)
        integer,intent(in)::types(:),mix(:,:)
        character(*),intent(in),optional::mode
        logical,intent(in),optional::diag
        real(dp),allocatable::g(:,:)
        character(len=12)::md
        logical::dg
        integer::n,nt,a,b,i,j,k,m,ne
        integer,allocatable::ri(:),rj(:),ord(:)
        n=size(types)
        nt=size(mix,1)
        md='digraph'
        if(present(mode))md=adjustl(mode)
        dg=.false.
        if(present(diag))dg=diag
        allocate(g(n,n))
        g=0
        if(trim(md)=='graph')then
            do a=1,nt
            do b=a,nt
            m=mix(a,b)
            ne=0
                do i=1,n
                do j=i,n
                if(types(i)/=a.or.types(j)/=b)cycle
                if(.not.dg.and.i==j)cycle
                if(a==b.and.types(j)/=a)cycle
                ne=ne+1
                end do
                end do
                if(m>ne)error stop 'rgnmix_exact: requested ties exceed available dyads'
                allocate(ri(max(1,ne)),rj(max(1,ne)),ord(max(1,ne)))
                k=0
                do i=1,n
                do j=i,n
                if(types(i)/=a.or.types(j)/=b)cycle
                if(.not.dg.and.i==j)cycle
                if(a==b.and.types(j)/=a)cycle
                k=k+1
                ri(k)=i
                rj(k)=j
                ord(k)=k
                end do
                end do
                if(ne>0)call shuffle_int(ord(:ne))
                do k=1,m
                i=ri(ord(k))
                j=rj(ord(k))
                g(i,j)=1
                g(j,i)=1
                end do
                deallocate(ri,rj,ord)
            end do
            end do
        else
            do a=1,nt
            do b=1,nt
            m=mix(a,b)
            ne=0
            do i=1,n
            do j=1,n
            if(types(i)==a.and.types(j)==b.and.(dg.or.i/=j))ne=ne+1
            end do
            end do
                if(m>ne)error stop 'rgnmix_exact: requested ties exceed available dyads'
                allocate(ri(max(1,ne)),rj(max(1,ne)),ord(max(1,ne)))
                k=0
                do i=1,n
                do j=1,n
                if(types(i)==a.and.types(j)==b.and.(dg.or.i/=j))then
                k=k+1
                ri(k)=i
                rj(k)=j
                ord(k)=k
                end if
                end do
                end do
                if(ne>0)call shuffle_int(ord(:ne))
                do k=1,m
                g(ri(ord(k)),rj(ord(k)))=1
                end do
                deallocate(ri,rj,ord)
            end do
            end do
        end if
    end function rgnmix_exact

    function rewire_ud(g,p) result(out)
        real(dp),intent(in)::g(:,:)
        real(dp),intent(in)::p
        real(dp),allocatable::out(:,:)
        integer::n,j,k,t,h
        real(dp)::u,tmp1,tmp2
        n=size(g,1)
        allocate(out(n,n))
        out=g
        if(n<3)return
        do j=1,n
        do k=j+1,n
        call random_number(u)
        if(u>=p)cycle
        t=j
        h=k
        call random_number(u)
            if(u<0.5_dp)then
            do
            h=randint(n)
            if(h/=j.and.h/=k)exit
            end do
            else
            do
            t=randint(n)
            if(t/=j.and.t/=k)exit
            end do
            end if
            tmp1=out(t,h)
            tmp2=out(h,t)
            out(t,h)=out(j,k)
            out(h,t)=out(k,j)
            out(j,k)=tmp1
            out(k,j)=tmp2
        end do
        end do
    end function rewire_ud

    function rewire_ws(g,p) result(out)
        real(dp),intent(in)::g(:,:)
        real(dp),intent(in)::p
        real(dp),allocatable::out(:,:)
        integer::n,j,k,t,h,tries
        real(dp)::u,tmp1,tmp2
        logical::ok
        n=size(g,1)
        allocate(out(n,n))
        out=g
        if(n<3)return
        do j=1,n
        do k=j+1,n
        if(g(j,k)==0.0_dp.and.g(k,j)==0.0_dp)cycle
        call random_number(u)
        if(u>=p)cycle
        ok=.false.
        tries=0
            do while(.not.ok.and.tries<10000)
            tries=tries+1
            t=j
            h=k
            call random_number(u)
            if(u<0.5_dp)then
            h=randint(n)
            else
            t=randint(n)
            end if
            ok=(h/=j.and.h/=k.and.t/=j.and.t/=k.and.out(t,h)==0.0_dp.and.out(h,t)==0.0_dp)
            end do
            if(.not.ok)cycle
            tmp1=out(t,h)
            tmp2=out(h,t)
            out(t,h)=out(j,k)
            out(h,t)=out(k,j)
            out(j,k)=tmp1
            out(k,j)=tmp2
        end do
        end do
    end function rewire_ws

    function rgws(nv,d,z,p) result(g)
        integer,intent(in)::nv,d,z
        real(dp),intent(in)::p
        real(dp),allocatable::g(:,:),lat(:,:)
        integer::tnv,i,j,dim,xi,xj,tmp,dist
        tnv=nv**d
        allocate(lat(tnv,tnv))
        lat=0
        do i=0,tnv-1
        do j=0,tnv-1
        if(i==j)cycle
        xi=i
        xj=j
        dist=0
        do dim=1,d
        dist=dist+abs(mod(xi,nv)-mod(xj,nv))
        xi=xi/nv
        xj=xj/nv
        end do
        if(dist<=z)lat(i+1,j+1)=1
        end do
        end do
        g=rewire_ws(lat,p)
    end function rgws

    subroutine shuffle_int(a)
        integer,intent(inout)::a(:)
        integer::i,j,t
        real(dp)::u
        do i=size(a),2,-1
        call random_number(u)
        j=1+int(u*real(i,dp))
        if(j>i)j=i
        t=a(i)
        a(i)=a(j)
        a(j)=t
        end do
    end subroutine shuffle_int

    integer function randint(n) result(k)
        integer,intent(in)::n
        real(dp)::u
        call random_number(u)
        k=1+int(u*real(n,dp))
        if(k>n)k=n
    end function randint


    function rgbn_mcmc(n_draws,nv,pi,sigma,rho,d,delta,epsilon,burn,thin,dichotomize_sib,max_density,seed_graph) result(draws)
        ! Biased-net MCMC corresponding to rgbn(method="mcmc")/bn_mcmc_R.
        integer,intent(in)::n_draws,nv
        real(dp),intent(in)::pi,sigma,rho,d(:,:)
        real(dp),intent(in),optional::delta,epsilon(:,:),max_density
        integer,intent(in),optional::burn,thin
        logical,intent(in),optional::dichotomize_sib
        real(dp),intent(in),optional::seed_graph(:,:)
        real(dp),allocatable::draws(:,:,:),g(:,:),eps(:,:)
        integer::b,th,total,step,save,j,k,x,parents,odeg
        real(dp)::del,md,u,ep,lnpar,lnsib,lndblr,lne,lni,lnsat
        logical::sd
        b=nv*nv*5*500
        if(present(burn))b=burn
        th=nv*nv*5
        if(present(thin))th=max(1,thin)
        del=0.0_dp
        if(present(delta))del=delta
        md=1.0_dp
        if(present(max_density))md=max_density
        sd=.false.
        if(present(dichotomize_sib))sd=dichotomize_sib
        allocate(g(nv,nv),eps(nv,nv),draws(n_draws,nv,nv))
        g=0.0_dp
        eps=0.0_dp
        draws=0.0_dp
        if(present(seed_graph))g=seed_graph
        if(present(epsilon))eps=epsilon
        do j=1,nv
        g(j,j)=0.0_dp
        end do
        lnpar=log(max(1.0_dp-pi,tiny(1.0_dp)))
        lnsib=log(max(1.0_dp-sigma,tiny(1.0_dp)))
        lndblr=log(max(1.0_dp-rho,tiny(1.0_dp)))
        lnsat=log(max(1.0_dp-del,tiny(1.0_dp)))
        total=b+n_draws*th
        save=0
        do step=1,total
            call random_number(u)
            j=1+int(u*real(nv,dp))
            if(j>nv)j=nv
            do
                call random_number(u)
                k=1+int(u*real(nv,dp))
                if(k>nv)k=nv
                if(k/=j)exit
            end do
            parents=0
            do x=1,nv
                if(x/=j.and.x/=k.and.g(x,j)>0.0_dp.and.g(x,k)>0.0_dp)parents=parents+1
            end do
            odeg=count(g(j,:)>0.0_dp)
            lne=log(max(1.0_dp-min(max(d(j,k),0.0_dp),1.0_dp),tiny(1.0_dp)))
            lni=log(max(1.0_dp-min(max(eps(j,k),0.0_dp),1.0_dp),tiny(1.0_dp)))
            if(sd)then
                ep=1.0_dp-exp(lne+merge(lnpar,0.0_dp,g(k,j)>0.0_dp)+merge(lnsib,0.0_dp,parents>0)+merge(lndblr,0.0_dp,g(k, &
                    & j)>0.0_dp.and.parents>0))
            else
                ep=1.0_dp-exp(lne+merge(lnpar,0.0_dp,g(k,j)>0.0_dp)+real(parents,dp)*lnsib+merge(real(parents,dp)*lndblr, &
                    & 0.0_dp,g(k,j)>0.0_dp))
            end if
            ep=ep*exp(real(odeg,dp)*lnsat+lni)
            ep=max(0.0_dp,min(1.0_dp,ep))
            call random_number(u)
            g(j,k)=merge(1.0_dp,0.0_dp,u<=ep)
            if(sum(g)>md*real(nv*(nv-1),dp))exit
            if(step>b.and.mod(step-b,th)==0)then
                save=save+1
                if(save<=n_draws)draws(save,:,:)=g
            end if
        end do
    end function rgbn_mcmc

end module sna_random
