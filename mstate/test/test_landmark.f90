program test_landmark
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use mstate
    implicit none
    type(transition_map)::tr
    type(msdata_type)::ms
    type(probtrans_type)::pt
    real(dp)::times(6,3),nanv
    integer::status(6,3),info,i
    call trans_illdeath(tr)
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    times(:,1)=nanv
    times(:,2)=[1.0_dp,1.0_dp,6.0_dp,6.0_dp,8.0_dp,9.0_dp]
    times(:,3)=[5.0_dp,1.0_dp,9.0_dp,7.0_dp,8.0_dp,12.0_dp]
    status(:,1)=0;status(:,2)=[1,0,1,1,0,1];status(:,3)=1
    call msprep(times,status,tr,ms,info=info)
    call landmark_aj(ms,tr,0.0_dp,[1],pt,info)
    call check(info==0 .and. pt%nt>1,'landmark result')
    do i=1,pt%nt
        call check(abs(sum(pt%p(i,1,:))-1.0_dp)<1e-10_dp,'landmark row sum')
    end do
    print '(a)','test_landmark: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
