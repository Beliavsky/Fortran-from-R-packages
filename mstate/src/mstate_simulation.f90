module mstate_simulation
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_is_finite
    use mstate_kinds, only : dp
    use mstate_types, only : transition_map, hazard_type, censor_distribution, simulated_msdata_type
    use mstate_transitions, only : enumerate_paths, is_circular
    implicit none
    private
    public :: simulate_state_probabilities, sample_path
    public :: mssample_state, mssample_paths, mssample_data, sample_path_general

contains

    subroutine sample_path(hz,tr,start_state,start_time,tvec,state_at,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::start_state
        real(dp),intent(in)::start_time,tvec(:)
        integer,intent(out)::state_at(:)
        integer,intent(out),optional::info
        call sample_path_general(hz,tr,start_state,start_time,tvec,state_at,info=info)
    end subroutine sample_path

    subroutine sample_path_general(hz,tr,start_state,start_time,tvec,state_at,clock,tstate,beta_state,cens,path_at,path_data,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::start_state
        real(dp),intent(in)::start_time,tvec(:)
        integer,intent(out)::state_at(:)
        character(len=*),intent(in),optional::clock
        real(dp),intent(in),optional::tstate(:),beta_state(:,:)
        type(censor_distribution),intent(in),optional::cens
        integer,intent(out),optional::path_at(:)
        type(simulated_msdata_type),intent(out),optional::path_data
        integer,intent(out),optional::info

        character(len=8)::clk
        real(dp),allocatable::ts(:)
        integer,allocatable::visited(:),paths(:,:),lens(:)
        integer::nvisit,cur,nr,chosen_trans,chosen_to,inf0
        real(dp)::tcur,tstop,base_cond,cens_time,cens_jump
        logical::censored,done
        type(simulated_msdata_type)::buf

        if(present(info))info=0
        if(size(state_at)/=size(tvec))then
            if(present(info))info=1;return
        end if
        if(start_state<1.or.start_state>tr%nstate)then
            if(present(info))info=2;return
        end if
        clk='forward';if(present(clock))clk=trim(adjustl(clock))
        if(clk/='forward'.and.clk/='reset')then
            if(present(info))info=3;return
        end if
        if(present(beta_state))then
            if(size(beta_state,1)/=tr%nstate.or.size(beta_state,2)/=tr%ntrans)then
                if(present(info))info=4;return
            end if
        end if
        if(present(tstate))then
            if(size(tstate)/=tr%nstate)then;if(present(info))info=5;return;end if
        end if
        if(present(path_at))then
            if(size(path_at)/=size(tvec))then;if(present(info))info=6;return;end if
            if(is_circular(tr))then;if(present(info))info=7;return;end if
        end if

        allocate(ts(tr%nstate));ts=0.0_dp
        if(present(tstate))ts=tstate
        allocate(visited(tr%nstate));visited=0
        visited(1)=start_state;nvisit=1
        state_at=0
        if(present(path_at))then
            call enumerate_paths(tr,start_state,paths,lens,inf0)
            if(inf0/=0)then;if(present(info))info=8;return;end if
            path_at=0
        end if
        call init_data_buffer(buf,max(16,2*max(1,tr%ntrans)))

        call sample_censor_once(cens,cens_time,cens_jump)
        cur=start_state;tcur=start_time;done=.false.
        do while(.not.done)
            nr=count(tr%trans(cur,:)>0)
            if(nr==0)then
                if(present(path_at))then
                    call fill_interval(tvec,state_at,path_at,paths,lens,visited,nvisit,tcur, &
                                       ieee_value(1.0_dp,ieee_positive_inf),cur)
                else
                    call fill_interval(tvec,state_at,t1=tcur,t2=ieee_value(1.0_dp,ieee_positive_inf),state=cur)
                end if
                exit
            end if
            if(clk=='forward')then;base_cond=tcur;else;base_cond=0.0_dp;end if
            call competing_risk_sample(hz,tr,cur,base_cond,ts,beta_state,cens_time,cens_jump, &
                                       tstop,chosen_trans,chosen_to,censored,inf0)
            if(inf0/=0.and.present(info))info=20+inf0
            if(clk=='reset'.and.ieee_is_finite(tstop))tstop=tstop+tcur

            if(present(path_at))then
                call fill_interval(tvec,state_at,path_at,paths,lens,visited,nvisit,tcur,tstop,cur)
            else
                call fill_interval(tvec,state_at,t1=tcur,t2=tstop,state=cur)
            end if
            call append_sojourn(buf,tr,cur,tcur,tstop,chosen_to)

            if(censored.or.chosen_trans==0.or..not.ieee_is_finite(tstop))then
                done=.true.
            else
                cur=chosen_to;tcur=tstop
                if(nvisit<tr%nstate)then
                    nvisit=nvisit+1;visited(nvisit)=cur
                end if
                if(present(beta_state))ts(cur)=tcur
            end if
        end do
        if(present(path_data))call finish_data_buffer(buf,path_data,1)
    end subroutine sample_path_general

    subroutine simulate_state_probabilities(hz,tr,m,tvec,pstate,start_state,start_time)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::m
        real(dp),intent(in)::tvec(:)
        real(dp),allocatable,intent(out)::pstate(:,:)
        integer,intent(in),optional::start_state
        real(dp),intent(in),optional::start_time
        integer::s0
        real(dp)::t0
        s0=1;if(present(start_state))s0=start_state
        t0=0.0_dp;if(present(start_time))t0=start_time
        call mssample_state(hz,tr,m,tvec,pstate,history_states=[s0],history_times=[t0])
    end subroutine simulate_state_probabilities

    subroutine mssample_state(hz,tr,m,tvec,pstate,history_states,history_times,history_tstates,beta_state,clock,cens,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::m
        real(dp),intent(in)::tvec(:)
        real(dp),allocatable,intent(out)::pstate(:,:)
        integer,intent(in),optional::history_states(:)
        real(dp),intent(in),optional::history_times(:),history_tstates(:,:),beta_state(:,:)
        character(len=*),intent(in),optional::clock
        type(censor_distribution),intent(in),optional::cens
        integer,intent(out),optional::info
        integer,allocatable::st(:)
        real(dp),allocatable::ts(:)
        integer::rep,i,s0,inf0
        real(dp)::t0
        if(present(info))info=0
        allocate(pstate(size(tvec),tr%nstate),st(size(tvec)),ts(tr%nstate));pstate=0.0_dp
        do rep=1,m
            call history_for_rep(rep,m,tr%nstate,history_states,history_times,history_tstates,s0,t0,ts,inf0)
            if(inf0/=0)then;if(present(info))info=inf0;return;end if
            call sample_path_general(hz,tr,s0,t0,tvec,st,clock=clock,tstate=ts,beta_state=beta_state,cens=cens,info=inf0)
            if(inf0/=0.and.present(info))info=inf0
            do i=1,size(tvec)
                if(st(i)>=1.and.st(i)<=tr%nstate)pstate(i,st(i))=pstate(i,st(i))+1.0_dp
            end do
        end do
        if(m>0)pstate=pstate/real(m,dp)
    end subroutine mssample_state

    subroutine mssample_paths(hz,tr,m,tvec,ppath,paths,lengths,history_states,history_times,history_tstates, &
                              beta_state,clock,cens,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::m
        real(dp),intent(in)::tvec(:)
        real(dp),allocatable,intent(out)::ppath(:,:)
        integer,allocatable,intent(out)::paths(:,:),lengths(:)
        integer,intent(in),optional::history_states(:)
        real(dp),intent(in),optional::history_times(:),history_tstates(:,:),beta_state(:,:)
        character(len=*),intent(in),optional::clock
        type(censor_distribution),intent(in),optional::cens
        integer,intent(out),optional::info
        integer,allocatable::st(:),pa(:)
        real(dp),allocatable::ts(:)
        integer::rep,i,s0,inf0,npath
        real(dp)::t0
        if(present(info))info=0
        s0=1;if(present(history_states))s0=history_states(1)
        call enumerate_paths(tr,s0,paths,lengths,inf0)
        if(inf0/=0)then;if(present(info))info=1;allocate(ppath(0,0));return;end if
        npath=size(lengths);allocate(ppath(size(tvec),npath),st(size(tvec)),pa(size(tvec)),ts(tr%nstate));ppath=0.0_dp
        do rep=1,m
            call history_for_rep(rep,m,tr%nstate,history_states,history_times,history_tstates,s0,t0,ts,inf0)
            if(inf0/=0)then;if(present(info))info=inf0;return;end if
            ! Upstream paths() is tied to the trajectory start state.  Mixed start states
            ! therefore need separate calls; reject them instead of silently mislabelling paths.
            if(s0/=paths(1,1))then;if(present(info))info=9;return;end if
            call sample_path_general(hz,tr,s0,t0,tvec,st,clock=clock,tstate=ts,beta_state=beta_state,cens=cens, &
                                     path_at=pa,info=inf0)
            do i=1,size(tvec)
                if(pa(i)>=1.and.pa(i)<=npath)ppath(i,pa(i))=ppath(i,pa(i))+1.0_dp
            end do
        end do
        if(m>0)ppath=ppath/real(m,dp)
    end subroutine mssample_paths

    subroutine mssample_data(hz,tr,m,data,history_states,history_times,history_tstates,beta_state,clock,cens,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::m
        type(simulated_msdata_type),intent(out)::data
        integer,intent(in),optional::history_states(:)
        real(dp),intent(in),optional::history_times(:),history_tstates(:,:),beta_state(:,:)
        character(len=*),intent(in),optional::clock
        type(censor_distribution),intent(in),optional::cens
        integer,intent(out),optional::info
        integer,allocatable::st(:)
        real(dp),allocatable::ts(:),dummy_t(:)
        type(simulated_msdata_type)::one,acc
        integer::rep,s0,inf0
        real(dp)::t0
        if(present(info))info=0
        allocate(dummy_t(1),st(1),ts(tr%nstate));dummy_t=0.0_dp
        call init_data_buffer(acc,max(16,m*max(1,tr%ntrans)))
        acc%n=0
        do rep=1,m
            call history_for_rep(rep,m,tr%nstate,history_states,history_times,history_tstates,s0,t0,ts,inf0)
            if(inf0/=0)then;if(present(info))info=inf0;return;end if
            call sample_path_general(hz,tr,s0,t0,dummy_t,st,clock=clock,tstate=ts,beta_state=beta_state,cens=cens, &
                                     path_data=one,info=inf0)
            call append_dataset(acc,one,rep)
        end do
        call trim_data(acc,data)
    end subroutine mssample_data

    subroutine competing_risk_sample(hz,tr,from_state,tcond,tstates,beta_state,cens_time,cens_jump, &
                                     event_time,transno,to_state,censored,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::from_state
        real(dp),intent(in)::tcond,tstates(:),cens_time,cens_jump
        real(dp),intent(out)::event_time
        integer,intent(out)::transno,to_state
        logical,intent(out)::censored
        real(dp),intent(in),optional::beta_state(:,:)
        integer,intent(out),optional::info
        integer,allocatable::trs(:),tos(:)
        real(dp),allocatable::scale(:),dh(:)
        integer::nr,j,i,idx,chosen
        real(dp)::u,sprev,snew,hjump,mass,cum,u2,eps_eq
        if(present(info))info=0
        nr=count(tr%trans(from_state,:)>0)
        allocate(trs(nr),tos(nr),scale(nr),dh(nr));j=0
        do i=1,tr%nstate
            if(tr%trans(from_state,i)>0)then
                j=j+1;trs(j)=tr%trans(from_state,i);tos(j)=i;scale(j)=1.0_dp
                if(present(beta_state))scale(j)=exp(sum(beta_state(:,trs(j))*tstates))
            end if
        end do
        call random_number(u);sprev=1.0_dp
        event_time=ieee_value(1.0_dp,ieee_positive_inf);transno=0;to_state=0;censored=.false.
        do idx=1,hz%nt
            if(hz%time(idx)<=tcond)cycle
            hjump=0.0_dp
            do j=1,nr
                if(idx==1)then
                    dh(j)=hz%haz(idx,trs(j))
                else
                    dh(j)=hz%haz(idx,trs(j))-hz%haz(idx-1,trs(j))
                end if
                dh(j)=max(0.0_dp,dh(j))*scale(j);hjump=hjump+dh(j)
            end do
            if(hjump<=0.0_dp)cycle
            if(hjump>1.0_dp)then
                if(present(info))info=1
                dh=dh/hjump;hjump=1.0_dp
            end if
            snew=sprev*(1.0_dp-hjump);mass=sprev-snew
            if(u> snew .and. u<=sprev)then
                event_time=hz%time(idx)
                if(cens_time<event_time)then
                    event_time=cens_time;censored=.true.;return
                end if
                eps_eq=32.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(cens_time),abs(event_time))
                if(abs(cens_time-event_time)<=eps_eq.and.ieee_is_finite(event_time))then
                    call random_number(u2);cum=0.0_dp;chosen=0
                    if(sum(dh)+max(0.0_dp,cens_jump)>0.0_dp)then
                        do j=1,nr
                            cum=cum+dh(j)/(sum(dh)+max(0.0_dp,cens_jump))
                            if(u2<=cum)then;chosen=j;exit;end if
                        end do
                    end if
                    if(chosen==0)then;censored=.true.;return;end if
                else
                    call random_number(u2);cum=0.0_dp;chosen=nr
                    do j=1,nr
                        cum=cum+dh(j)/sum(dh)
                        if(u2<=cum)then;chosen=j;exit;end if
                    end do
                end if
                transno=trs(chosen);to_state=tos(chosen);return
            end if
            sprev=snew
        end do
        if(cens_time<event_time)then;event_time=cens_time;censored=.true.;end if
    end subroutine competing_risk_sample

    subroutine sample_censor_once(cens,time,jump)
        type(censor_distribution),intent(in),optional::cens
        real(dp),intent(out)::time,jump
        real(dp),allocatable::prob(:)
        real(dp)::u,cum,fprev,fcur
        integer::i,n
        time=ieee_value(1.0_dp,ieee_positive_inf);jump=0.0_dp
        if(.not.present(cens))return
        n=cens%n;if(n<=0)return
        allocate(prob(n+1));prob=0.0_dp;fprev=0.0_dp
        do i=1,n
            if(allocated(cens%surv))then
                if(size(cens%surv)>=n)then
                    fcur=1.0_dp-cens%surv(i)
                else
                    fcur=fprev
                end if
            else if(allocated(cens%haz))then
                if(size(cens%haz)>=n)then
                    fcur=1.0_dp-exp(-cens%haz(i))
                else
                    fcur=fprev
                end if
            else
                fcur=fprev
            end if
            prob(i)=max(0.0_dp,fcur-fprev);fprev=fcur
        end do
        prob(n+1)=max(0.0_dp,1.0_dp-fprev)
        if(sum(prob)<=0.0_dp)return
        prob=prob/sum(prob);call random_number(u);cum=0.0_dp
        do i=1,n+1
            cum=cum+prob(i)
            if(u<=cum)exit
        end do
        if(i<=n)then
            time=cens%time(i)
            if(allocated(cens%haz))then
                if(size(cens%haz)>=n)then
                    if(i==1)then;jump=cens%haz(1);else;jump=cens%haz(i)-cens%haz(i-1);end if
                end if
            end if
        end if
    end subroutine sample_censor_once

    subroutine fill_interval(tvec,state_at,path_at,paths,lens,visited,nvisit,t1,t2,state)
        real(dp),intent(in)::tvec(:),t1,t2
        integer,intent(inout)::state_at(:)
        integer,intent(in),optional::paths(:,:),lens(:),visited(:),nvisit
        integer,intent(inout),optional::path_at(:)
        integer,intent(in)::state
        integer::i,pidx
        pidx=0
        if(present(path_at).and.present(paths).and.present(lens).and.present(visited).and.present(nvisit))then
            pidx=find_prefix(paths,lens,visited,nvisit)
        end if
        do i=1,size(tvec)
            if(tvec(i)>=t1.and.tvec(i)<t2)then
                state_at(i)=state
                if(present(path_at).and.pidx>0)path_at(i)=pidx
            end if
        end do
    end subroutine fill_interval

    integer function find_prefix(paths,lens,visited,nvisit) result(idx)
        integer,intent(in)::paths(:,:),lens(:),visited(:),nvisit
        integer::i
        idx=0
        do i=1,size(lens)
            if(lens(i)==nvisit)then
                if(all(paths(i,1:nvisit)==visited(1:nvisit)))then;idx=i;return;end if
            end if
        end do
    end function find_prefix

    subroutine history_for_rep(rep,m,k,states,times,tstates,state,time,ts,info)
        integer,intent(in)::rep,m,k
        integer,intent(in),optional::states(:)
        real(dp),intent(in),optional::times(:),tstates(:,:)
        integer,intent(out)::state,info
        real(dp),intent(out)::time,ts(k)
        info=0;state=1;time=0.0_dp;ts=0.0_dp
        if(present(states))then
            if(size(states)==1)then;state=states(1)
            else if(size(states)==m)then;state=states(rep)
            else;info=10;return;end if
        end if
        if(present(times))then
            if(size(times)==1)then;time=times(1)
            else if(size(times)==m)then;time=times(rep)
            else;info=11;return;end if
        end if
        if(present(tstates))then
            if(size(tstates,1)/=k)then;info=12;return;end if
            if(size(tstates,2)==1)then;ts=tstates(:,1)
            else if(size(tstates,2)==m)then;ts=tstates(:,rep)
            else;info=13;return;end if
        end if
    end subroutine history_for_rep

    subroutine init_data_buffer(x,cap)
        type(simulated_msdata_type),intent(out)::x
        integer,intent(in)::cap
        integer::c
        c=max(1,cap);x%n=0
        allocate(x%id(c),x%from(c),x%to(c),x%status(c),x%trans(c),x%tstart(c),x%tstop(c),x%duration(c))
        x%id=0;x%from=0;x%to=0;x%status=0;x%trans=0;x%tstart=0.0_dp;x%tstop=0.0_dp;x%duration=0.0_dp
    end subroutine init_data_buffer

    subroutine ensure_capacity(x,need)
        type(simulated_msdata_type),intent(inout)::x
        integer,intent(in)::need
        integer::old,newc
        old=size(x%id);if(need<=old)return
        newc=max(need,2*old)
        call grow_int(x%id,newc);call grow_int(x%from,newc);call grow_int(x%to,newc)
        call grow_int(x%status,newc);call grow_int(x%trans,newc)
        call grow_real(x%tstart,newc);call grow_real(x%tstop,newc);call grow_real(x%duration,newc)
    end subroutine ensure_capacity

    subroutine grow_int(a,nnew)
        integer,allocatable,intent(inout)::a(:)
        integer,intent(in)::nnew
        integer,allocatable::b(:)
        allocate(b(nnew));b=0;b(1:size(a))=a;call move_alloc(b,a)
    end subroutine grow_int

    subroutine grow_real(a,nnew)
        real(dp),allocatable,intent(inout)::a(:)
        integer,intent(in)::nnew
        real(dp),allocatable::b(:)
        allocate(b(nnew));b=0.0_dp;b(1:size(a))=a;call move_alloc(b,a)
    end subroutine grow_real

    subroutine append_sojourn(buf,tr,from,tstart,tstop,chosen_to)
        type(simulated_msdata_type),intent(inout)::buf
        type(transition_map),intent(in)::tr
        integer,intent(in)::from,chosen_to
        real(dp),intent(in)::tstart,tstop
        integer::j,nr,pos
        nr=count(tr%trans(from,:)>0);if(nr==0)return
        call ensure_capacity(buf,buf%n+nr)
        do j=1,tr%nstate
            if(tr%trans(from,j)<=0)cycle
            buf%n=buf%n+1;pos=buf%n
            buf%from(pos)=from;buf%to(pos)=j;buf%trans(pos)=tr%trans(from,j)
            buf%status(pos)=merge(1,0,j==chosen_to);buf%tstart(pos)=tstart;buf%tstop(pos)=tstop
            buf%duration(pos)=tstop-tstart
        end do
    end subroutine append_sojourn

    subroutine finish_data_buffer(buf,out,idval)
        type(simulated_msdata_type),intent(in)::buf
        type(simulated_msdata_type),intent(out)::out
        integer,intent(in)::idval
        integer::n
        n=buf%n;out%n=n
        allocate(out%id(n),out%from(n),out%to(n),out%status(n),out%trans(n),out%tstart(n),out%tstop(n),out%duration(n))
        if(n>0)then
            out%id=idval;out%from=buf%from(1:n);out%to=buf%to(1:n);out%status=buf%status(1:n);out%trans=buf%trans(1:n)
            out%tstart=buf%tstart(1:n);out%tstop=buf%tstop(1:n);out%duration=buf%duration(1:n)
        end if
    end subroutine finish_data_buffer

    subroutine append_dataset(acc,one,idval)
        type(simulated_msdata_type),intent(inout)::acc
        type(simulated_msdata_type),intent(in)::one
        integer,intent(in)::idval
        integer::a,b,n
        n=one%n;if(n==0)return
        a=acc%n+1;b=acc%n+n;call ensure_capacity(acc,b)
        acc%id(a:b)=idval;acc%from(a:b)=one%from;acc%to(a:b)=one%to;acc%status(a:b)=one%status;acc%trans(a:b)=one%trans
        acc%tstart(a:b)=one%tstart;acc%tstop(a:b)=one%tstop;acc%duration(a:b)=one%duration;acc%n=b
    end subroutine append_dataset

    subroutine trim_data(acc,out)
        type(simulated_msdata_type),intent(in)::acc
        type(simulated_msdata_type),intent(out)::out
        call finish_data_buffer(acc,out,0)
        if(out%n>0)out%id=acc%id(1:acc%n)
    end subroutine trim_data


end module mstate_simulation
