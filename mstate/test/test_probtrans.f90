program test_probtrans
    use mstate
    implicit none
    type(transition_map)::tr
    type(hazard_type)::hz
    type(probtrans_type)::pt
    real(dp),allocatable::elos(:,:)
    integer::info
    call trans_comprisk(2,tr)
    hz%nt=2;hz%ntrans=2
    allocate(hz%time(2),hz%haz(2,2),hz%varhaz(2,2,2))
    hz%time=[1.0_dp,2.0_dp]
    hz%haz=reshape([0.2_dp,0.3_dp,0.1_dp,0.3_dp],[2,2])
    hz%varhaz=0.0_dp
    hz%varhaz(1,1,1)=0.02_dp;hz%varhaz(1,2,2)=0.01_dp
    hz%varhaz(2,1,1)=0.03_dp;hz%varhaz(2,2,2)=0.03_dp
    call probtrans(hz,tr,0.0_dp,pt,direction='forward',method='greenwood',variance=.true.,info=info)
    call check(info==0 .and. pt%nt==3,'probtrans dimensions')
    call close(pt%p(2,1,1),0.7_dp,1e-12_dp,'p11 t1')
    call close(pt%p(2,1,2),0.2_dp,1e-12_dp,'p12 t1')
    call close(pt%p(3,1,1),0.49_dp,1e-12_dp,'p11 t2')
    call close(pt%p(3,1,2),0.27_dp,1e-12_dp,'p12 t2')
    call close(pt%p(3,1,3),0.24_dp,1e-12_dp,'p13 t2')
    call check(all(pt%se>=0.0_dp),'nonnegative se')
    call expected_length_of_stay(pt,2.0_dp,elos,info)
    call close(elos(1,1),1.7_dp,1e-12_dp,'elos 11')
    call close(elos(1,2),0.2_dp,1e-12_dp,'elos 12')
    call close(elos(1,3),0.1_dp,1e-12_dp,'elos 13')
    call close(elos(2,2),2.0_dp,1e-12_dp,'absorbing elos')
    print '(a)','test_probtrans: PASS'
contains
    subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if;end subroutine
    subroutine close(x,y,tol,msg);real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
        call check(abs(x-y)<=tol,msg);end subroutine
end program
