program test_nonparametric_sim
    use mstate
    implicit none
    real(dp)::time(6)
    integer::status(6),info
    real(dp),allocatable::ut(:),surv(:),cif(:,:),pstate(:,:)
    type(transition_map)::tr
    type(hazard_type)::hz
    real(dp)::tvec(3)
    time=[1.0_dp,2.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    status=[1,2,0,1,0,2]
    call cumulative_incidence(time,status,ut,surv,cif,info)
    call check(info==0 .and. size(ut)==4,'cuminc event times')
    call check(all(cif>=0.0_dp),'cuminc nonnegative')
    call check(all(surv>=0.0_dp .and. surv<=1.0_dp),'survival range')
    call trans_comprisk(2,tr)
    hz%nt=2;hz%ntrans=2;allocate(hz%time(2),hz%haz(2,2),hz%varhaz(2,2,2))
    hz%time=[1.0_dp,2.0_dp];hz%haz=reshape([0.2_dp,0.3_dp,0.1_dp,0.3_dp],[2,2]);hz%varhaz=0.0_dp
    tvec=[0.5_dp,1.5_dp,2.5_dp]
    call simulate_state_probabilities(hz,tr,3000,tvec,pstate)
    call check(all(abs(sum(pstate,dim=2)-1.0_dp)<1e-12_dp),'sim probabilities sum')
    call check(pstate(1,1)>0.99_dp,'before first hazard')
    print '(a)','test_nonparametric_sim: PASS'
contains
    subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if;end subroutine
end program
