program test_msfit
    use mstate
    implicit none
    real(dp)::start(3),stop(3),score(3),times(2),newrisk(1)
    integer::event(3),strata(2),kstrata(1),info
    real(dp),allocatable::xmat(:,:),vcov(:,:),newx(:,:)
    type(hazard_type)::hz, hz2
    real(dp), allocatable :: beta(:)
    start=0.0_dp;stop=[1.0_dp,2.0_dp,3.0_dp];event=[1,1,0];score=1.0_dp
    strata=[1,4];kstrata=1;times=[1.0_dp,2.0_dp];newrisk=1.0_dp
    allocate(xmat(3,0),vcov(0,0),newx(1,0))
    call agmssurv(start,stop,event,score,xmat,vcov,strata,kstrata,times,newx,newrisk,1,.true.,hz,info)
    call check(info==0,'agmssurv info')
    call close(hz%haz(1,1),1.0_dp/3.0_dp,1e-12_dp,'breslow h1')
    call close(hz%haz(2,1),5.0_dp/6.0_dp,1e-12_dp,'breslow h2')
    call close(hz%varhaz(2,1,1),13.0_dp/36.0_dp,1e-12_dp,'breslow variance')
    allocate(beta(0))
    call msfit_from_cox_arrays(start,stop,event,xmat,beta,vcov,strata,kstrata,times,newx,1,.true.,hz2,info)
    call close(hz2%haz(2,1),hz%haz(2,1),1e-12_dp,'cox-array wrapper')
    print '(a)','test_msfit: PASS'
contains
    subroutine check(ok,msg);logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if;end subroutine
    subroutine close(x,y,tol,msg);real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
        call check(abs(x-y)<=tol,msg);end subroutine
end program
