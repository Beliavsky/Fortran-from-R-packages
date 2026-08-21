program test_survival_bridge
    use mstate
    use survival_cox, only : coxph_fit_counting
    implicit none
    real(dp) :: start1(5), stop1(5), x1(5,1)
    real(dp) :: start2(10), stop2(10), x2(10,1), times(3), newx(2,1)
    integer :: status1(5), status2(10), strata(10), kstrata(2), info
    type(coxph_result) :: one, two, fit
    type(hazard_type) :: hz

    start1 = 0.0_dp
    stop1 = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    status1 = [1,1,0,1,0]
    x1(:,1) = [-1.0_dp,0.5_dp,1.0_dp,-0.5_dp,0.2_dp]
    call coxph_fit_counting(start1,stop1,status1,x1,one,'breslow')

    start2(1:5)=start1;start2(6:10)=start1
    stop2(1:5)=stop1;stop2(6:10)=stop1
    status2(1:5)=status1;status2(6:10)=status1
    x2(1:5,:)=x1;x2(6:10,:)=x1
    strata=[1,1,1,1,1,2,2,2,2,2]
    call coxph_fit_stratified_counting(start2,stop2,status2,x2,strata,two,'breslow')
    call close(two%coef(1),one%coef(1),1.0e-9_dp,'duplicated strata coefficient')
    call close(two%var(1,1),0.5_dp*one%var(1,1),1.0e-8_dp,'duplicated strata variance')

    times=[1.0_dp,2.0_dp,4.0_dp];newx=0.0_dp;kstrata=[1,2]
    call msfit_cox(start2,stop2,status2,x2,strata,newx,kstrata,times,1,.true.,fit,hz,info)
    call check(info==0,'msfit_cox info')
    call check(hz%ntrans==2.and.hz%nt==3,'msfit_cox shape')
    call close(hz%haz(3,1),hz%haz(3,2),1.0e-10_dp,'identical stratum hazards')
    print '(a)','test_survival_bridge: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
    subroutine close(x,y,tol,msg)
        real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
        call check(abs(x-y)<=tol,msg)
    end subroutine
end program test_survival_bridge
