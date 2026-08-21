program test_relative
    use mstate
    implicit none
    type(transition_map)::tr,trn
    type(hazard_type)::hz,hzn
    integer,allocatable::link(:,:),nlink(:),smap(:,:)
    integer::info
    real(dp)::pop(1,2)
    call trans_comprisk(2,tr)
    call modify_transition_relative(tr,[2],trn,link,nlink,smap,info)
    call check(info==0.and.trn%nstate==4.and.trn%ntrans==3,'relative transition dimensions')
    call check(nlink(1)==1.and.nlink(2)==2,'relative transition linkage')
    hz%nt=1;hz%ntrans=2
    allocate(hz%time(1),hz%haz(1,2),hz%varhaz(1,2,2))
    hz%time=1.0_dp;hz%haz(1,:)=[0.2_dp,0.3_dp];hz%varhaz=0.0_dp
    hz%varhaz(1,1,1)=0.01_dp;hz%varhaz(1,2,2)=0.04_dp
    hz%varhaz(1,1,2)=0.005_dp;hz%varhaz(1,2,1)=0.005_dp
    pop=0.0_dp;pop(1,2)=0.1_dp
    call split_relative_hazards(hz,tr,[2],pop,hzn,trn,link,nlink,info)
    call check(info==0,'relative hazard split')
    call check(abs(hzn%haz(1,link(1,2))-0.1_dp)<1e-12_dp,'population hazard')
    call check(abs(hzn%haz(1,link(2,2))-0.2_dp)<1e-12_dp,'excess hazard')
    call check(abs(hzn%varhaz(1,link(1,2),link(1,2)))<1e-12_dp,'population variance fixed')
    call check(abs(hzn%varhaz(1,link(2,2),link(2,2))-0.04_dp)<1e-12_dp,'excess variance')
    call check(abs(hzn%varhaz(1,link(1,1),link(2,2))-0.005_dp)<1e-12_dp,'excess covariance')
    print '(a)','test_relative: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
