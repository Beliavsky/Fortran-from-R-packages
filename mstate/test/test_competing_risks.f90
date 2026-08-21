program test_competing_risks
    use mstate
    implicit none
    real(dp)::time(4)
    integer::status(4),group(4),info
    type(cuminc_result)::fit
    type(cuminc_result),allocatable::fits(:)
    integer,allocatable::groups(:)
    time=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    status=[1,2,1,2]
    call cumulative_incidence_fit(time,status,fit,info)
    call check(info==0 .and. fit%nt==4 .and. fit%ncause==2,'fit dimensions')
    call check(abs(fit%surv(4))<1.0e-12_dp,'final survival')
    call check(maxval(abs(fit%cif(4,:)-[0.5_dp,0.5_dp]))<1.0e-12_dp,'final CIF')
    call check(maxval(abs(fit%se_cif(4,:)-0.25_dp))<1.0e-10_dp,'binomial final SE')
    call check(all(fit%se_cif>=0.0_dp).and.all(fit%se_surv>=0.0_dp),'nonnegative SE')
    group=[1,1,2,2]
    call cumulative_incidence_grouped(time,status,group,groups,fits,info)
    call check(info==0.and.size(groups)==2.and.all(groups==[1,2]),'group labels')
    call check(size(fits)==2.and.all([(fits(info)%nt==2,info=1,2)]),'group curves')
    print '(a)','test_competing_risks: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
