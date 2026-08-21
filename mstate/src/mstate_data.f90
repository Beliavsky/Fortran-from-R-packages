module mstate_data
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use mstate_kinds, only : dp
    use mstate_types, only : transition_map, msdata_type
    use mstate_transitions, only : is_circular
    implicit none
    private
    public :: msprep, cut_landmark, xsect, event_counts, expand_covariates

contains

    subroutine msprep(times, status, tr, ms, ids, start_state, start_time, info)
        real(dp), intent(in) :: times(:, :)
        integer, intent(in) :: status(:, :)
        type(transition_map), intent(in) :: tr
        type(msdata_type), intent(out) :: ms
        integer, intent(in), optional :: ids(:), start_state(:)
        real(dp), intent(in), optional :: start_time(:)
        integer, intent(out), optional :: info
        integer :: n, s, i, j, cur, nxt, nr, k, row, maxrows
        integer, allocatable :: tos(:), trs(:), aid(:), afrom(:), ato(:), atr(:), astat(:)
        real(dp), allocatable :: ats(:), ate(:)
        real(dp) :: t0, tevent, tcens, tnext
        logical :: have_event

        if (present(info)) info=0
        n=size(times,1); s=size(times,2)
        if (size(status,1)/=n .or. size(status,2)/=s .or. s/=tr%nstate) then
            if (present(info)) info=1
            return
        end if
        if (is_circular(tr)) then
            if (present(info)) info=2
            return
        end if
        if (present(ids)) then
            if (size(ids)/=n) then; if(present(info))info=3; return; end if
        end if
        if (present(start_state)) then
            if (size(start_state)/=n) then; if(present(info))info=4; return; end if
        end if
        if (present(start_time)) then
            if (size(start_time)/=n) then; if(present(info))info=5; return; end if
        end if

        maxrows=max(1,n*max(1,tr%ntrans))
        allocate(aid(maxrows),afrom(maxrows),ato(maxrows),atr(maxrows),astat(maxrows))
        allocate(ats(maxrows),ate(maxrows)); row=0

        do i=1,n
            cur=1; if(present(start_state)) cur=start_state(i)
            t0=0.0_dp; if(present(start_time)) t0=start_time(i)
            do k=1,tr%nstate+1
                nr=count(tr%trans(cur,:)>0)
                if(nr==0) exit
                allocate(tos(nr),trs(nr)); nr=0
                do j=1,tr%nstate
                    if(tr%trans(cur,j)>0) then
                        nr=nr+1; tos(nr)=j; trs(nr)=tr%trans(cur,j)
                    end if
                end do
                tevent=huge(1.0_dp); tcens=huge(1.0_dp); nxt=0; have_event=.false.
                do j=1,nr
                    if(ieee_is_finite(times(i,tos(j))) .and. times(i,tos(j))>=t0) then
                        tcens=min(tcens,times(i,tos(j)))
                        if(status(i,tos(j))==1 .and. times(i,tos(j))<tevent) then
                            tevent=times(i,tos(j)); nxt=tos(j); have_event=.true.
                        end if
                    end if
                end do
                if(have_event) then
                    tnext=tevent
                else
                    tnext=tcens
                end if
                if(.not.ieee_is_finite(tnext) .or. tnext>=huge(1.0_dp)/2) then
                    deallocate(tos,trs); exit
                end if
                do j=1,nr
                    row=row+1
                    aid(row)=i; if(present(ids)) aid(row)=ids(i)
                    afrom(row)=cur; ato(row)=tos(j); atr(row)=trs(j)
                    ats(row)=t0; ate(row)=tnext
                    astat(row)=merge(1,0,have_event .and. tos(j)==nxt)
                end do
                deallocate(tos,trs)
                if(.not.have_event) exit
                cur=nxt; t0=tnext
            end do
        end do
        ms%n=row
        allocate(ms%id(row),ms%from(row),ms%to(row),ms%trans(row),ms%status(row))
        allocate(ms%tstart(row),ms%tstop(row),ms%time(row))
        if(row>0) then
            ms%id=aid(1:row); ms%from=afrom(1:row); ms%to=ato(1:row)
            ms%trans=atr(1:row); ms%status=astat(1:row)
            ms%tstart=ats(1:row); ms%tstop=ate(1:row); ms%time=ate(1:row)-ats(1:row)
            call sort_msdata(ms)
        end if
    end subroutine msprep

    subroutine sort_msdata(ms)
        type(msdata_type),intent(inout)::ms
        integer::i,j,ki,kf,kt,ktr,ks
        real(dp)::ka,kb,km
        do i=2,ms%n
            ki=ms%id(i); kf=ms%from(i); kt=ms%to(i); ktr=ms%trans(i); ks=ms%status(i)
            ka=ms%tstart(i); kb=ms%tstop(i); km=ms%time(i); j=i-1
            do while(j>=1)
                if(ms%id(j)<ki) exit
                if(ms%id(j)==ki .and. ms%tstart(j)<ka) exit
                if(ms%id(j)==ki .and. ms%tstart(j)==ka .and. ms%from(j)<kf) exit
                if(ms%id(j)==ki .and. ms%tstart(j)==ka .and. ms%from(j)==kf .and. ms%to(j)<=kt) exit
                ms%id(j+1)=ms%id(j); ms%from(j+1)=ms%from(j); ms%to(j+1)=ms%to(j)
                ms%trans(j+1)=ms%trans(j); ms%status(j+1)=ms%status(j)
                ms%tstart(j+1)=ms%tstart(j); ms%tstop(j+1)=ms%tstop(j); ms%time(j+1)=ms%time(j); j=j-1
            end do
            ms%id(j+1)=ki; ms%from(j+1)=kf; ms%to(j+1)=kt; ms%trans(j+1)=ktr; ms%status(j+1)=ks
            ms%tstart(j+1)=ka; ms%tstop(j+1)=kb; ms%time(j+1)=km
        end do
    end subroutine sort_msdata

    subroutine cut_landmark(ms, lm, out, cens)
        type(msdata_type),intent(in)::ms
        real(dp),intent(in)::lm
        type(msdata_type),intent(out)::out
        real(dp),intent(in),optional::cens
        integer::i,n
        n=0
        do i=1,ms%n
            if(ms%tstop(i)>lm) then
                if(.not.present(cens) .or. min(ms%tstop(i),cens)>=max(ms%tstart(i),lm)) n=n+1
            end if
        end do
        out%n=n
        allocate(out%id(n),out%from(n),out%to(n),out%trans(n),out%status(n),out%tstart(n),out%tstop(n),out%time(n))
        n=0
        do i=1,ms%n
            if(ms%tstop(i)<=lm) cycle
            if(present(cens)) then
                if(min(ms%tstop(i),cens)<max(ms%tstart(i),lm)) cycle
            end if
            n=n+1; out%id(n)=ms%id(i); out%from(n)=ms%from(i); out%to(n)=ms%to(i); out%trans(n)=ms%trans(i)
            out%tstart(n)=max(ms%tstart(i),lm); out%tstop(n)=ms%tstop(i); out%status(n)=ms%status(i)
            if(present(cens)) then
                if(out%tstop(n)>cens) then; out%tstop(n)=cens; out%status(n)=0; end if
            end if
            out%time(n)=out%tstop(n)-out%tstart(n)
        end do
    end subroutine cut_landmark

    subroutine xsect(ms, xtime, ids, states)
        type(msdata_type),intent(in)::ms
        real(dp),intent(in)::xtime
        integer,allocatable,intent(out)::ids(:),states(:)
        integer,allocatable::ti(:),ts(:)
        integer::i,j,n
        allocate(ti(ms%n),ts(ms%n)); n=0
        do i=1,ms%n
            if(ms%tstart(i)<=xtime .and. ms%tstop(i)>xtime) then
                j=0
                if(n>0) then
                    if(any(ti(1:n)==ms%id(i))) j=1
                end if
                if(j==0) then; n=n+1; ti(n)=ms%id(i); ts(n)=ms%from(i); end if
            end if
        end do
        allocate(ids(n),states(n)); if(n>0) then; ids=ti(1:n); states=ts(1:n); end if
    end subroutine xsect

    subroutine event_counts(ms,tr,counts,total_entering)
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        integer,allocatable,intent(out)::counts(:,:),total_entering(:)
        integer::i,s
        allocate(counts(tr%nstate,tr%nstate),total_entering(tr%nstate)); counts=0; total_entering=0
        do i=1,ms%n
            if(ms%status(i)==1) counts(ms%from(i),ms%to(i))=counts(ms%from(i),ms%to(i))+1
        end do
        do s=1,tr%nstate
            do i=1,ms%n
                if(ms%from(i)==s) total_entering(s)=max(total_entering(s),count([(ms%id==ms%id(i)).and.(ms%from==s)]))
            end do
            ! For nonabsorbing states, each entrant contributes one row per outgoing transition.
            if(count(tr%trans(s,:)>0)>0) then
                total_entering(s)=count([(ms%from(i)==s, i=1,ms%n)])/count(tr%trans(s,:)>0)
            else
                total_entering(s)=sum(counts(:,s))
            end if
        end do
    end subroutine event_counts

    subroutine expand_covariates(cov, transno, ntrans, expanded)
        real(dp),intent(in)::cov(:,:)
        integer,intent(in)::transno(:),ntrans
        real(dp),allocatable,intent(out)::expanded(:,:)
        integer::n,p,i,j,k
        n=size(cov,1); p=size(cov,2)
        if(size(transno)/=n) then; allocate(expanded(0,0)); return; end if
        allocate(expanded(n,p*ntrans)); expanded=0.0_dp
        do i=1,n
            k=transno(i)
            if(k<1.or.k>ntrans) cycle
            do j=1,p
                expanded(i,(j-1)*ntrans+k)=cov(i,j)
            end do
        end do
    end subroutine expand_covariates
end module mstate_data
