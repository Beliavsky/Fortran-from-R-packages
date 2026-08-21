program test_efron
    use mstate
    implicit none
    real(dp)::start(3),stop(3),score(3),times(1),newrisk(1)
    integer::event(3),strata(2),kstrata(1),info
    real(dp),allocatable::xmat(:,:),vcov(:,:),newx(:,:)
    type(hazard_type)::hb,he
    start=0.0_dp;stop=1.0_dp;event=[1,1,0];score=1.0_dp
    strata=[1,4];kstrata=1;times=1.0_dp;newrisk=1.0_dp
    allocate(xmat(3,0),vcov(0,0),newx(1,0))
    call agmssurv(start,stop,event,score,xmat,vcov,strata,kstrata,times,newx,newrisk,1,.true.,hb,info)
    call check(info==0,'breslow info')
    call agmssurv(start,stop,event,score,xmat,vcov,strata,kstrata,times,newx,newrisk,2,.true.,he,info)
    call check(info==0,'efron info')
    call close(hb%haz(1,1),2.0_dp/3.0_dp,1e-12_dp,'breslow tied hazard')
    call close(he%haz(1,1),5.0_dp/6.0_dp,1e-12_dp,'efron tied hazard')
    call close(he%varhaz(1,1,1),13.0_dp/36.0_dp,1e-12_dp,'efron tied variance')
    print '(a)','test_efron: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
    subroutine close(x,y,tol,msg)
        real(dp),intent(in)::x,y,tol;character(len=*),intent(in)::msg
        call check(abs(x-y)<=tol,msg)
    end subroutine
end program
