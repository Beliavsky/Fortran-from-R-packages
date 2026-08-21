program test_utilities
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use mstate
    implicit none
    type(transition_map)::tr
    type(msdata_type)::ms,ms2,boot
    type(etmdata_type)::etm
    real(dp)::times(6,3),nanv,grid(2)
    real(dp),allocatable::w(:),wm(:,:)
    integer::status(6,3),info
    call trans_illdeath(tr)
    nanv=ieee_value(0.0_dp,ieee_quiet_nan)
    times(:,1)=nanv
    times(:,2)=[1.0_dp,1.0_dp,6.0_dp,6.0_dp,8.0_dp,9.0_dp]
    times(:,3)=[5.0_dp,1.0_dp,9.0_dp,7.0_dp,8.0_dp,12.0_dp]
    status(:,1)=0;status(:,2)=[1,0,1,1,0,1];status(:,3)=1
    call msprep(times,status,tr,ms,info=info)
    call msdata_to_etm(ms,etm)
    call check(etm%n==10,'etm interval count')
    call etm_to_msdata(etm,tr,ms2,info)
    call check(info==0 .and. ms2%n==ms%n,'etm roundtrip rows')
    call check(count(ms2%status==1)==count(ms%status==1),'etm roundtrip events')
    call bootstrap_msdata(ms,boot)
    call check(boot%n>0,'bootstrap nonempty')
    grid=[0.5_dp,4.0_dp]
    call optimal_weights_multiple(ms,tr,grid,1,w,0.0_dp,info)
    call check(info==0 .and. size(w)==2 .and. all(w>=0.0_dp),'optimal multiple weights')
    call optimal_weights_matrix(ms,tr,grid,1,wm,0.0_dp,info)
    call check(info==0 .and. all(shape(wm)==[2,3]) .and. all(wm>=0.0_dp),'optimal matrix weights')
    print '(a)','test_utilities: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
