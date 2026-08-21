program test_mssample_modes
    use mstate
    implicit none
    type(transition_map)::tr,tr2
    type(hazard_type)::hz,hz2
    type(censor_distribution)::cens
    type(simulated_msdata_type)::dat
    real(dp)::tvec(2),beta(2,1),tstates(2,1)
    real(dp),allocatable::pstate(:,:),ppath(:,:)
    integer,allocatable::st(:),paths(:,:),lens(:)
    integer::info

    call trans_comprisk(1,tr)
    hz%nt=1;hz%ntrans=1
    allocate(hz%time(1),hz%haz(1,1),hz%varhaz(1,1,1))
    hz%time=2.0_dp;hz%haz=1.0_dp;hz%varhaz=0.0_dp
    tvec=[5.5_dp,7.5_dp];allocate(st(2))
    call sample_path_general(hz,tr,1,5.0_dp,tvec,st,clock='forward',info=info)
    call check(info==0.and.all(st==1),'forward clock after hazard support')
    call sample_path_general(hz,tr,1,5.0_dp,tvec,st,clock='reset',path_data=dat,info=info)
    call check(info==0.and.all(st==[1,2]),'reset clock transition')
    call check(dat%n==1.and.dat%status(1)==1.and.abs(dat%tstop(1)-7.0_dp)<1e-12_dp,'reset data output')

    ! History-state effects multiply the relevant baseline cumulative hazard.
    hz%time=1.0_dp;hz%haz=0.2_dp
    beta=0.0_dp;beta(1,1)=log(5.0_dp);tstates(:,1)=[1.0_dp,0.0_dp]
    tvec=[0.5_dp,1.5_dp]
    call mssample_state(hz,tr,1,tvec,pstate,history_tstates=tstates,beta_state=beta,info=info)
    call check(info==0.and.pstate(1,1)>0.999_dp.and.pstate(2,2)>0.999_dp,'history-state hazard effect')

    ! Independent censoring is sampled once and terminates the current sojourn.
    cens%n=1;allocate(cens%time(1),cens%haz(1),cens%surv(1))
    cens%time=0.5_dp;cens%haz=1.0_dp;cens%surv=0.0_dp
    hz%time=1.0_dp;hz%haz=1.0_dp;tvec=[0.25_dp,0.75_dp]
    call sample_path_general(hz,tr,1,0.0_dp,tvec,st,cens=cens,path_data=dat,info=info)
    call check(info==0.and.st(1)==1.and.st(2)==0,'censoring state output')
    call check(dat%n==1.and.dat%status(1)==0.and.abs(dat%tstop(1)-0.5_dp)<1e-12_dp,'censoring data output')

    ! Path and long-data outputs follow upstream paths(), including prefixes.
    call trans_comprisk(2,tr2)
    hz2%nt=1;hz2%ntrans=2
    allocate(hz2%time(1),hz2%haz(1,2),hz2%varhaz(1,2,2))
    hz2%time=1.0_dp;hz2%haz(1,:)=[1.0_dp,0.0_dp];hz2%varhaz=0.0_dp
    tvec=[0.5_dp,1.5_dp]
    call mssample_paths(hz2,tr2,4,tvec,ppath,paths,lens,info=info)
    call check(info==0.and.size(lens)==3.and.all(lens==[1,2,2]),'path prefix catalogue')
    call check(ppath(1,1)>0.999_dp.and.ppath(2,2)>0.999_dp,'path probabilities')
    call mssample_data(hz2,tr2,1,dat,info=info)
    call check(info==0.and.dat%n==2.and.sum(dat%status)==1,'long data competing rows')
    call check(all(dat%to==[2,3]).and.all(dat%trans==[1,2]),'long data destinations')
    print '(a)','test_mssample_modes: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
