module mstate_relative
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use mstate_kinds, only : dp
    use mstate_types, only : transition_map, hazard_type, msdata_type, relative_bootstrap_type, relative_msfit_type
    use mstate_transitions, only : transition_from_matrix
    use mstate_nonparametric, only : nelson_aalen_msdata
    use relsurv_ratetable, only : ratetable_type, net_summary_type, expprep2_summary
    implicit none
    private
    public :: modify_transition_relative, split_relative_hazards
    public :: haz_function_relsurv, population_hazards_relsurv
    public :: msfit_relsurv, msboot_relsurv, msfit_relsurv_full, hazard_add_times

contains

    subroutine modify_transition_relative(tr,split_trans,tr_new,link,nlink,state_map,info)
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(transition_map),intent(out)::tr_new
        integer,allocatable,intent(out)::link(:,:),nlink(:),state_map(:,:)
        integer,intent(out),optional::info
        integer,allocatable::mode(:),base(:),pop(:),exc(:),mat(:,:)
        integer::j,r,nnew,src,dst,newsrc,newdst,val,nsplit,nincoming,inf0
        logical::issplit
        if(present(info))info=0
        allocate(mode(tr%nstate),base(tr%nstate),pop(tr%nstate),exc(tr%nstate));mode=0;base=0;pop=0;exc=0
        do j=1,tr%nstate
            nincoming=count(tr%to==j)
            nsplit=0
            do r=1,tr%ntrans
                if(tr%to(r)==j.and.any(split_trans==r))nsplit=nsplit+1
            end do
            if(nsplit>0.and.any(tr%trans(j,:)>0))then
                if(present(info))info=1;allocate(link(0,0),nlink(0),state_map(0,0));return
            end if
            if(nsplit==0)then
                mode(j)=1
            else if(nsplit==nincoming)then
                mode(j)=2
            else
                mode(j)=3
            end if
        end do
        nnew=0
        do j=1,tr%nstate
            select case(mode(j))
            case(1)
                nnew=nnew+1;base(j)=nnew
            case(2)
                nnew=nnew+1;pop(j)=nnew
                nnew=nnew+1;exc(j)=nnew
            case(3)
                nnew=nnew+1;base(j)=nnew
                nnew=nnew+1;pop(j)=nnew
                nnew=nnew+1;exc(j)=nnew
            end select
        end do
        allocate(mat(nnew,nnew));mat=0
        allocate(link(2,tr%ntrans),nlink(tr%ntrans),state_map(3,tr%nstate));link=0;nlink=0
        state_map(1,:)=base;state_map(2,:)=pop;state_map(3,:)=exc
        val=0
        ! R's modify_transMat renumbers transitions in row-major order. Build all
        ! destination edges first, then assign the same row-major numbering.
        do r=1,tr%ntrans
            src=tr%from(r);dst=tr%to(r)
            newsrc=base(src)
            if(newsrc==0)then
                if(present(info))info=2;return
            end if
            issplit=any(split_trans==r)
            if(issplit)then
                if(pop(dst)==0.or.exc(dst)==0)then;if(present(info))info=3;return;end if
                mat(newsrc,pop(dst))=-r
                mat(newsrc,exc(dst))=-r
            else
                if(mode(dst)==3)then;newdst=base(dst)
                else if(mode(dst)==1)then;newdst=base(dst)
                else
                    if(present(info))info=4;return
                end if
                mat(newsrc,newdst)=-r
            end if
        end do
        do src=1,nnew
            do dst=1,nnew
                if(mat(src,dst)/=0)then
                    r=-mat(src,dst);val=val+1;mat(src,dst)=val
                    nlink(r)=nlink(r)+1;link(nlink(r),r)=val
                end if
            end do
        end do
        call transition_from_matrix(mat,tr_new,inf0)
        if(inf0/=0.and.present(info))info=10+inf0
    end subroutine modify_transition_relative

    subroutine split_relative_hazards(hz,tr,split_trans,pop_haz,hz_new,tr_new,link,nlink,info)
        type(hazard_type),intent(in)::hz
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        real(dp),intent(in)::pop_haz(:,:)
        type(hazard_type),intent(out)::hz_new
        type(transition_map),intent(out)::tr_new
        integer,allocatable,intent(out)::link(:,:),nlink(:)
        integer,intent(out),optional::info
        integer,allocatable::state_map(:,:),orig_of(:),kind(:)
        integer::r,q,a,b,inf0
        if(present(info))info=0
        if(hz%ntrans/=tr%ntrans.or.size(pop_haz,1)/=hz%nt.or.size(pop_haz,2)/=tr%ntrans)then
            if(present(info))info=1;allocate(link(0,0),nlink(0));return
        end if
        call modify_transition_relative(tr,split_trans,tr_new,link,nlink,state_map,inf0)
        if(inf0/=0)then;if(present(info))info=2;return;end if
        hz_new%nt=hz%nt;hz_new%ntrans=tr_new%ntrans
        allocate(hz_new%time(hz%nt),hz_new%haz(hz%nt,tr_new%ntrans))
        allocate(hz_new%varhaz(hz%nt,tr_new%ntrans,tr_new%ntrans))
        hz_new%time=hz%time;hz_new%haz=0.0_dp;hz_new%varhaz=0.0_dp
        allocate(orig_of(tr_new%ntrans),kind(tr_new%ntrans));orig_of=0;kind=0
        do r=1,tr%ntrans
            if(nlink(r)==1)then
                q=link(1,r);hz_new%haz(:,q)=hz%haz(:,r);orig_of(q)=r;kind(q)=0
            else if(nlink(r)==2)then
                q=link(1,r);hz_new%haz(:,q)=pop_haz(:,r);orig_of(q)=r;kind(q)=1
                q=link(2,r);hz_new%haz(:,q)=hz%haz(:,r)-pop_haz(:,r);orig_of(q)=r;kind(q)=2
            end if
        end do
        if(allocated(hz%varhaz))then
            do a=1,tr_new%ntrans
                do b=1,tr_new%ntrans
                    if(kind(a)==1.or.kind(b)==1)cycle
                    hz_new%varhaz(:,a,b)=hz%varhaz(:,orig_of(a),orig_of(b))
                end do
            end do
        end if
    end subroutine split_relative_hazards

    subroutine haz_function_relsurv(ms,xrate,transition,tab,add_times,eval_times,haz_pop,haz_excess, &
                                    nrisk,nevent,ncensor,se,info,time_scale,allow_no_events)
        type(msdata_type),intent(in)::ms
        real(dp),intent(in)::xrate(:,:)
        integer,intent(in)::transition
        type(ratetable_type),intent(in)::tab
        real(dp),intent(in)::add_times(:)
        real(dp),allocatable,intent(out)::eval_times(:),haz_pop(:),haz_excess(:),nrisk(:),nevent(:),ncensor(:),se(:)
        integer,intent(out),optional::info
        real(dp),intent(in),optional::time_scale
        logical,intent(in),optional::allow_no_events
        type(net_summary_type)::summ
        real(dp),allocatable::x(:,:),y(:),ys(:),grid(:)
        integer,allocatable::stat(:)
        real(dp)::scale,den,acc
        integer::i,j,n
        logical::allow0

        if(present(info))info=0
        if(size(xrate,1)/=ms%n.or.size(xrate,2)/=tab%ndim)then
            if(present(info))info=1
            call empty_haz_function(eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
            return
        end if
        if(transition<1)then
            if(present(info))info=2
            call empty_haz_function(eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
            return
        end if
        n=count(ms%trans==transition)
        allow0=.false.;if(present(allow_no_events))allow0=allow_no_events
        if(n==0)then
            if(allow0)then
                call copy_zero_grid(add_times,eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
                return
            end if
            if(present(info))info=3
            call empty_haz_function(eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
            return
        end if
        if(.not.allow0.and.count((ms%trans==transition).and.(ms%status==1))==0)then
            if(present(info))info=4
            call empty_haz_function(eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
            return
        end if
        scale=1.0_dp;if(present(time_scale))scale=time_scale
        if(scale<=0.0_dp)then
            if(present(info))info=5
            call empty_haz_function(eval_times,haz_pop,haz_excess,nrisk,nevent,ncensor,se)
            return
        end if
        allocate(x(n,tab%ndim),y(n),ys(n),stat(n));j=0
        do i=1,ms%n
            if(ms%trans(i)/=transition)cycle
            j=j+1;x(j,:)=xrate(i,:);y(j)=ms%tstop(i)*scale;ys(j)=ms%tstart(i)*scale;stat(j)=ms%status(i)
        end do
        call merged_time_grid(add_times*scale,y,grid)
        call expprep2_summary(tab,x,y,stat,grid,summ,fast=.true.,ystart=ys)
        allocate(eval_times(size(grid)),haz_pop(size(grid)),haz_excess(size(grid)),nrisk(size(grid)), &
                 nevent(size(grid)),ncensor(size(grid)),se(size(grid)))
        eval_times=grid/scale;haz_pop=0.0_dp;haz_excess=0.0_dp;nrisk=summ%yi;nevent=summ%dni;ncensor=0.0_dp;se=0.0_dp
        acc=0.0_dp
        do j=1,size(grid)
            den=summ%yi(j)
            if(den>0.0_dp.and.ieee_is_finite(den))then
                haz_pop(j)=summ%yidli(j)/den
                haz_excess(j)=summ%dni(j)/den-haz_pop(j)
                acc=acc+summ%dni(j)/(den*den)
            end if
            se(j)=sqrt(max(0.0_dp,acc))
            if(j<size(grid))then
                ncensor(j)=-summ%yi(j+1)+summ%yi(j)-summ%dni(j)
            else
                ncensor(j)=summ%yi(j)-summ%dni(j)
            end if
        end do
    end subroutine haz_function_relsurv

    subroutine empty_haz_function(t,hp,he,nr,ne,nc,se)
        real(dp),allocatable,intent(out)::t(:),hp(:),he(:),nr(:),ne(:),nc(:),se(:)
        allocate(t(0),hp(0),he(0),nr(0),ne(0),nc(0),se(0))
    end subroutine empty_haz_function

    subroutine copy_zero_grid(add_times,t,hp,he,nr,ne,nc,se)
        real(dp),intent(in)::add_times(:)
        real(dp),allocatable,intent(out)::t(:),hp(:),he(:),nr(:),ne(:),nc(:),se(:)
        allocate(t(size(add_times)),hp(size(add_times)),he(size(add_times)),nr(size(add_times)), &
                 ne(size(add_times)),nc(size(add_times)),se(size(add_times)))
        t=add_times;hp=0.0_dp;he=0.0_dp;nr=0.0_dp;ne=0.0_dp;nc=0.0_dp;se=0.0_dp
    end subroutine copy_zero_grid

    subroutine population_hazards_relsurv(ms,tr,split_trans,tab,xrate,times,pop_haz,info,time_scale,allow_no_events)
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(ratetable_type),intent(in)::tab
        real(dp),intent(in)::xrate(:,:),times(:)
        real(dp),allocatable,intent(out)::pop_haz(:,:)
        integer,intent(out),optional::info
        real(dp),intent(in),optional::time_scale
        logical,intent(in),optional::allow_no_events
        real(dp),allocatable::hp(:),he(:),nr(:),ne(:),nc(:),se(:)
        integer::r,ierr
        logical::allow0
        if(present(info))info=0
        allocate(pop_haz(size(times),tr%ntrans));pop_haz=0.0_dp
        allow0=.false.;if(present(allow_no_events))allow0=allow_no_events
        do r=1,tr%ntrans
            if(.not.any(split_trans==r))cycle
            block
                real(dp),allocatable::gt(:),cum(:)
                integer::i,j
                call haz_function_relsurv(ms,xrate,r,tab,times,gt,hp,he,nr,ne,nc,se,ierr,time_scale,allow0)
                if(ierr/=0)then
                    if(present(info))info=100*r+ierr
                    return
                end if
                allocate(cum(size(gt)));cum=0.0_dp
                do i=1,size(gt)
                    cum(i)=hp(i)
                    if(i>1)cum(i)=cum(i)+cum(i-1)
                end do
                do i=1,size(times)
                    j=find_time(gt,times(i))
                    if(j>0)pop_haz(i,r)=cum(j)
                end do
            end block
        end do
    end subroutine population_hazards_relsurv

    subroutine msfit_relsurv(hz,ms,tr,split_trans,tab,xrate,hz_new,tr_new,link,nlink,info,time_scale, &
                             allow_no_events,pop_haz_out)
        type(hazard_type),intent(in)::hz
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(ratetable_type),intent(in)::tab
        real(dp),intent(in)::xrate(:,:)
        type(hazard_type),intent(out)::hz_new
        type(transition_map),intent(out)::tr_new
        integer,allocatable,intent(out)::link(:,:),nlink(:)
        integer,intent(out),optional::info
        real(dp),intent(in),optional::time_scale
        logical,intent(in),optional::allow_no_events
        real(dp),allocatable,intent(out),optional::pop_haz_out(:,:)
        real(dp),allocatable::pop(:,:)
        integer::ierr
        logical::allow0
        if(present(info))info=0
        allow0=.false.;if(present(allow_no_events))allow0=allow_no_events
        if(hz%ntrans/=tr%ntrans.or.size(xrate,1)/=ms%n)then;if(present(info))info=1;return;end if
        call population_hazards_relsurv(ms,tr,split_trans,tab,xrate,hz%time,pop,ierr,time_scale,allow0)
        if(ierr/=0)then;if(present(info))info=ierr;return;end if
        call split_relative_hazards(hz,tr,split_trans,pop,hz_new,tr_new,link,nlink,ierr)
        if(ierr/=0)then;if(present(info))info=10000+ierr;return;end if
        if(present(pop_haz_out))then
            allocate(pop_haz_out(size(pop,1),size(pop,2)));pop_haz_out=pop
        end if
    end subroutine msfit_relsurv

    subroutine msboot_relsurv(ms,tr,split_trans,tab,xrate,times,b,result,info,time_scale,seed,boot_original)
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(ratetable_type),intent(in)::tab
        real(dp),intent(in)::xrate(:,:),times(:)
        integer,intent(in)::b
        type(relative_bootstrap_type),intent(out)::result
        integer,intent(out),optional::info
        real(dp),intent(in),optional::time_scale
        integer,intent(in),optional::seed
        logical,intent(in),optional::boot_original
        type(msdata_type)::bms
        type(hazard_type)::bhz,bhz_grid,brs
        type(transition_map)::tr_new
        real(dp),allocatable::bx(:,:),vals(:)
        integer,allocatable::link(:,:),nlink(:),src_ids(:)
        logical,allocatable::valid(:,:),valid_orig(:,:)
        logical::keep_original
        integer::ib,r,q,ierr,nv
        real(dp)::nanv

        if(present(info))info=0
        if(b<2.or.size(xrate,1)/=ms%n)then;if(present(info))info=1;return;end if
        if(present(seed))call set_seed_integer(seed)
        call modify_transition_relative_only(tr,split_trans,tr_new,link,nlink,ierr)
        if(ierr/=0)then;if(present(info))info=2;return;end if
        result%nt=size(times);result%ntrans=tr_new%ntrans;result%norig=tr%ntrans;result%b=b
        keep_original=.false.;if(present(boot_original))keep_original=boot_original
        allocate(result%haz(result%nt,result%ntrans,b),result%varhaz(result%nt,result%ntrans), &
                 result%nvalid(result%ntrans),result%valid_rep(b),valid(result%ntrans,b))
        nanv=ieee_value(0.0_dp,ieee_quiet_nan);result%haz=nanv;result%varhaz=nanv;result%nvalid=0
        result%valid_rep=.false.;valid=.false.;result%has_original=keep_original
        if(keep_original)then
            allocate(result%original_haz(result%nt,tr%ntrans,b),result%original_varhaz(result%nt,tr%ntrans), &
                     result%original_nvalid(tr%ntrans),valid_orig(tr%ntrans,b))
            result%original_haz=nanv;result%original_varhaz=nanv;result%original_nvalid=0;valid_orig=.false.
        end if
        do ib=1,b
            call bootstrap_msdata_rate(ms,xrate,bms,bx,src_ids)
            if(size(src_ids)==0.or.count_distinct(src_ids)<2)cycle
            result%valid_rep(ib)=.true.
            call nelson_aalen_msdata(bms,tr,bhz,info=ierr)
            if(ierr/=0)cycle
            call project_hazard_grid(bhz,times,bhz_grid)
            if(keep_original)then
                do r=1,tr%ntrans
                    if(.not.any(bms%trans==r))cycle
                    valid_orig(r,ib)=.true.;result%original_haz(:,r,ib)=bhz_grid%haz(:,r)
                end do
            end if
            call msfit_relsurv(bhz_grid,bms,tr,split_trans,tab,bx,brs,tr_new,link,nlink,ierr,time_scale, &
                               allow_no_events=.true.)
            if(ierr/=0)cycle
            do r=1,tr%ntrans
                if(.not.any(bms%trans==r))cycle
                do q=1,nlink(r)
                    valid(link(q,r),ib)=.true.
                    result%haz(:,link(q,r),ib)=brs%haz(:,link(q,r))
                end do
            end do
        end do
        do q=1,result%ntrans
            result%nvalid(q)=count(valid(q,:))
            if(result%nvalid(q)<2)cycle
            allocate(vals(result%nvalid(q)))
            do r=1,result%nt
                nv=0
                do ib=1,b
                    if(valid(q,ib))then;nv=nv+1;vals(nv)=result%haz(r,q,ib);end if
                end do
                result%varhaz(r,q)=sample_variance(vals)
            end do
            deallocate(vals)
        end do
        if(keep_original)then
            do q=1,tr%ntrans
                result%original_nvalid(q)=count(valid_orig(q,:))
                if(result%original_nvalid(q)<2)cycle
                allocate(vals(result%original_nvalid(q)))
                do r=1,result%nt
                    nv=0
                    do ib=1,b
                        if(valid_orig(q,ib))then;nv=nv+1;vals(nv)=result%original_haz(r,q,ib);end if
                    end do
                    result%original_varhaz(r,q)=sample_variance(vals)
                end do
                deallocate(vals)
            end do
        end if
    end subroutine msboot_relsurv

    subroutine hazard_add_times(hz,add_times,out)
        type(hazard_type),intent(in)::hz
        real(dp),intent(in)::add_times(:)
        type(hazard_type),intent(out)::out
        real(dp),allocatable::grid(:)
        integer::i,j,k
        call merged_time_grid(hz%time,add_times,grid)
        out%nt=size(grid);out%ntrans=hz%ntrans
        allocate(out%time(out%nt),out%haz(out%nt,out%ntrans),out%varhaz(out%nt,out%ntrans,out%ntrans))
        out%time=grid;out%haz=0.0_dp;out%varhaz=0.0_dp
        do i=1,out%nt
            j=0
            do k=1,hz%nt
                if(hz%time(k)<=grid(i)+time_tol(hz%time(k),grid(i)))j=k
            end do
            if(j>0)then
                out%haz(i,:)=hz%haz(j,:)
                if(allocated(hz%varhaz))out%varhaz(i,:,:)=hz%varhaz(j,:,:)
            end if
        end do
    end subroutine hazard_add_times

    subroutine msfit_relsurv_full(hz,ms,tr,split_trans,tab,xrate,result,variance_mode,b,seed,time_format,info,add_times)
        type(hazard_type),intent(in)::hz
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(ratetable_type),intent(in)::tab
        real(dp),intent(in)::xrate(:,:)
        type(relative_msfit_type),intent(out)::result
        character(len=*),intent(in),optional::variance_mode,time_format
        integer,intent(in),optional::b,seed
        integer,intent(out),optional::info
        real(dp),intent(in),optional::add_times(:)
        character(len=9)::mode
        real(dp)::scale
        integer::nb,ierr,q
        type(hazard_type)::workhz

        if(present(info))info=0
        mode='fixed';if(present(variance_mode))mode=trim(adjustl(variance_mode))
        if(mode/='fixed'.and.mode/='bootstrap'.and.mode/='both')then;if(present(info))info=1;return;end if
        scale=1.0_dp
        if(present(time_format))then
            select case(trim(adjustl(time_format)))
            case('days');scale=1.0_dp
            case('years');scale=365.241_dp
            case('months');scale=365.241_dp/12.0_dp
            case default;if(present(info))info=2;return
            end select
        end if
        if(present(add_times))then
            call hazard_add_times(hz,add_times,workhz)
        else
            workhz=hz
        end if
        call msfit_relsurv(workhz,ms,tr,split_trans,tab,xrate,result%fit,result%trans,result%link,result%nlink, &
                           ierr,time_scale=scale,pop_haz_out=result%population_haz)
        if(ierr/=0)then;if(present(info))info=100+ierr;return;end if
        result%variance_mode=mode;result%has_bootstrap=.false.
        if(mode=='bootstrap'.or.mode=='both')then
            nb=10;if(present(b))nb=b
            if(present(seed))then
                call msboot_relsurv(ms,tr,split_trans,tab,xrate,workhz%time,nb,result%bootstrap,ierr, &
                                    time_scale=scale,seed=seed)
            else
                call msboot_relsurv(ms,tr,split_trans,tab,xrate,workhz%time,nb,result%bootstrap,ierr,time_scale=scale)
            end if
            if(ierr/=0)then;if(present(info))info=200+ierr;return;end if
            if(any(result%bootstrap%nvalid<2))then;if(present(info))info=3;return;end if
            result%has_bootstrap=.true.;result%bootstrap_fit=result%fit;result%bootstrap_fit%varhaz=0.0_dp
            do q=1,result%bootstrap_fit%ntrans
                result%bootstrap_fit%varhaz(:,q,q)=result%bootstrap%varhaz(:,q)
            end do
        end if
    end subroutine msfit_relsurv_full

    subroutine modify_transition_relative_only(tr,split_trans,tr_new,link,nlink,info)
        type(transition_map),intent(in)::tr
        integer,intent(in)::split_trans(:)
        type(transition_map),intent(out)::tr_new
        integer,allocatable,intent(out)::link(:,:),nlink(:)
        integer,intent(out)::info
        integer,allocatable::state_map(:,:)
        call modify_transition_relative(tr,split_trans,tr_new,link,nlink,state_map,info)
    end subroutine modify_transition_relative_only

    subroutine project_hazard_grid(src,times,dst)
        type(hazard_type),intent(in)::src
        real(dp),intent(in)::times(:)
        type(hazard_type),intent(out)::dst
        integer::i,j,k
        dst%nt=size(times);dst%ntrans=src%ntrans
        allocate(dst%time(dst%nt),dst%haz(dst%nt,dst%ntrans),dst%varhaz(dst%nt,dst%ntrans,dst%ntrans))
        dst%time=times;dst%haz=0.0_dp;dst%varhaz=0.0_dp
        do i=1,dst%nt
            j=0
            do k=1,src%nt
                if(src%time(k)<=times(i)+time_tol(src%time(k),times(i)))j=k
            end do
            if(j>0)then;dst%haz(i,:)=src%haz(j,:);dst%varhaz(i,:,:)=src%varhaz(j,:,:);end if
        end do
    end subroutine project_hazard_grid

    subroutine bootstrap_msdata_rate(ms,xrate,out,xout,source_ids)
        type(msdata_type),intent(in)::ms
        real(dp),intent(in)::xrate(:,:)
        type(msdata_type),intent(out)::out
        real(dp),allocatable,intent(out)::xout(:,:)
        integer,allocatable,intent(out)::source_ids(:)
        integer,allocatable::uid(:),pick(:),cnt(:)
        integer::i,j,k,nuid,nrow,row
        real(dp)::u
        call unique_ids_local(ms%id,uid);nuid=size(uid)
        allocate(pick(nuid),cnt(nuid),source_ids(nuid));cnt=0
        do i=1,nuid
            call random_number(u);k=min(nuid,1+int(u*real(nuid,dp)));pick(i)=uid(k);source_ids(i)=pick(i)
            cnt(i)=count(ms%id==pick(i))
        end do
        nrow=sum(cnt);out%n=nrow
        allocate(out%id(nrow),out%from(nrow),out%to(nrow),out%trans(nrow),out%status(nrow))
        allocate(out%tstart(nrow),out%tstop(nrow),out%time(nrow),xout(nrow,size(xrate,2)))
        row=0
        do i=1,nuid
            do j=1,ms%n
                if(ms%id(j)/=pick(i))cycle
                row=row+1;out%id(row)=i;out%from(row)=ms%from(j);out%to(row)=ms%to(j);out%trans(row)=ms%trans(j)
                out%status(row)=ms%status(j);out%tstart(row)=ms%tstart(j);out%tstop(row)=ms%tstop(j);out%time(row)=ms%time(j)
                xout(row,:)=xrate(j,:)
            end do
        end do
    end subroutine bootstrap_msdata_rate

    subroutine merged_time_grid(a,b,g)
        real(dp),intent(in)::a(:),b(:)
        real(dp),allocatable,intent(out)::g(:)
        real(dp),allocatable::tmp(:)
        integer::i,n
        allocate(tmp(size(a)+size(b)));tmp=[a,b];call sort_real_local(tmp)
        if(size(tmp)==0)then;allocate(g(0));return;end if
        n=1
        do i=2,size(tmp)
            if(.not.same_time_local(tmp(i),tmp(n)))then;n=n+1;tmp(n)=tmp(i);end if
        end do
        allocate(g(n));g=tmp(1:n)
    end subroutine merged_time_grid

    subroutine sort_real_local(x)
        real(dp),intent(inout)::x(:)
        integer::i,j
        real(dp)::key
        do i=2,size(x)
            key=x(i);j=i-1
            do while(j>=1)
                if(x(j)<=key)exit
                x(j+1)=x(j);j=j-1
            end do
            x(j+1)=key
        end do
    end subroutine sort_real_local

    integer function find_time(x,v) result(idx)
        real(dp),intent(in)::x(:),v
        integer::i
        idx=0
        do i=1,size(x)
            if(same_time_local(x(i),v))then;idx=i;return;end if
        end do
    end function find_time

    pure logical function same_time_local(a,b)
        real(dp),intent(in)::a,b
        same_time_local=abs(a-b)<=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(b))
    end function same_time_local

    pure real(dp) function time_tol(a,b)
        real(dp),intent(in)::a,b
        time_tol=64.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(a),abs(b))
    end function time_tol

    subroutine unique_ids_local(ids,u)
        integer,intent(in)::ids(:)
        integer,allocatable,intent(out)::u(:)
        integer,allocatable::tmp(:)
        integer::i,n
        allocate(tmp(size(ids)));n=0
        do i=1,size(ids)
            if(n==0)then;n=1;tmp(1)=ids(i)
            else if(.not.any(tmp(1:n)==ids(i)))then;n=n+1;tmp(n)=ids(i)
            end if
        end do
        allocate(u(n));if(n>0)u=tmp(1:n)
    end subroutine unique_ids_local

    integer function count_distinct(x) result(n)
        integer,intent(in)::x(:)
        integer,allocatable::u(:)
        call unique_ids_local(x,u);n=size(u)
    end function count_distinct

    pure real(dp) function sample_variance(x) result(v)
        real(dp),intent(in)::x(:)
        real(dp)::m
        if(size(x)<2)then;v=0.0_dp;return;end if
        m=sum(x)/real(size(x),dp);v=sum((x-m)**2)/real(size(x)-1,dp)
    end function sample_variance

    subroutine set_seed_integer(seed)
        integer,intent(in)::seed
        integer,allocatable::put(:)
        integer::i,n
        call random_seed(size=n);allocate(put(n))
        do i=1,n;put(i)=mod(abs(seed)+104729*i,2147483646)+1;end do
        call random_seed(put=put)
    end subroutine set_seed_integer

end module mstate_relative
