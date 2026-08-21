program test_transitions_data
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use mstate
    implicit none
    type(transition_map) :: tr, tr2
    type(msdata_type) :: ms, cut
    real(dp) :: times(6,3)
    integer :: status(6,3), info
    integer, allocatable :: q(:,:), absv(:), paths(:,:), lens(:), ids(:), states(:), counts(:,:), total(:)
    logical, allocatable :: tra(:,:)
    real(dp), allocatable :: cov(:,:), ex(:,:)
    real(dp) :: nanv
    integer :: i

    call trans_illdeath(tr)
    call check(tr%nstate==3 .and. tr%ntrans==3,'illdeath dimensions')
    call check(all(tr%trans==reshape([0,0,0,1,0,0,2,3,0],[3,3])),'illdeath transition numbers')
    call trans2q(tr,q)
    call check(q(1,3)>0,'transitive reachability')
    call absorbing_states(tr,absv); call check(size(absv)==1 .and. absv(1)==3,'absorbing')
    call check(.not.is_circular(tr),'acyclic')
    call trans2tra(tr,tra); call tra2trans(tra,tr2,info)
    call check(info==0.and.all(tr2%trans==tr%trans),'etm transition round trip')
    call enumerate_paths(tr,1,paths,lens,info)
    call check(info==0 .and. size(lens)==4,'path prefix count')

    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    times(:,1)=nanv
    times(:,2)=[1.0_dp,1.0_dp,6.0_dp,6.0_dp,8.0_dp,9.0_dp]
    times(:,3)=[5.0_dp,1.0_dp,9.0_dp,7.0_dp,8.0_dp,12.0_dp]
    status(:,1)=0
    status(:,2)=[1,0,1,1,0,1]
    status(:,3)=1
    call msprep(times,status,tr,ms,info=info)
    call check(info==0,'msprep info')
    call check(ms%n==16,'msprep row count')
    call check(count(ms%status==1)==10,'msprep event count')
    call check(count(ms%trans==3)==4,'post-illness risk rows')
    call event_counts(ms,tr,counts,total)
    call check(counts(1,2)==4 .and. counts(1,3)==2 .and. counts(2,3)==4,'event table')
    call check(all(total==[6,4,6]),'entering counts')
    call xsect(ms,0.5_dp,ids,states)
    call check(size(ids)==6 .and. all(states==1),'xsect at 0.5')
    call cut_landmark(ms,6.5_dp,cut,cens=10.0_dp)
    call check(all(cut%tstart>=6.5_dp) .and. all(cut%tstop<=10.0_dp),'landmark cut')

    allocate(cov(ms%n,2)); cov(:,1)=real(ms%id,dp);cov(:,2)=1.0_dp
    call expand_covariates(cov,ms%trans,tr%ntrans,ex)
    call check(size(ex,2)==6,'expanded cov dimension')
    call check(all([(count(abs(ex(i,:))>0.0_dp)==2,i=1,ms%n)]),'expanded cov sparsity')
    print '(a)','test_transitions_data: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
