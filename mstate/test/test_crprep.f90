program test_crprep
    use mstate
    implicit none
    real(dp)::tstop(6),tstart(6),keep(6,1)
    integer::status(6),strata(6),ids(6),info,i
    type(crprep_type)::w,wt
    tstop=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
    status=[1,2,0,1,2,1]
    strata=[1,1,1,1,2,2]
    ids=[11,12,13,14,15,16]
    keep(:,1)=[10.0_dp,20.0_dp,30.0_dp,40.0_dp,50.0_dp,60.0_dp]
    call crprep(tstop,status,[1],0,w,id=ids,strata=strata,keep=keep,shorten=.false.,info=info)
    call check(info==0.and.w%n==8.and.w%nkeep==1,'expanded row count')
    call check(count(w%id==12)==2.and.count(w%id==15)==2,'competing-event extension')
    call check(maxval(abs(pack(w%tstop,w%id==12)-[2.0_dp,4.0_dp]))<1e-12_dp,'subject 12 stop grid')
    call check(maxval(abs(pack(w%weight_cens,w%id==12)-[1.0_dp,0.5_dp]))<1e-12_dp,'censor weights')
    call check(all(pack(w%weight_cens,w%strata==2)==1.0_dp),'stratum-specific censoring')
    call check(all(pack(w%keep(:,1),w%id==15)==50.0_dp),'kept covariate replication')
    tstart=[0.0_dp,0.5_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp]
    call crprep(tstop,status,[1],0,wt,tstart=tstart,id=ids,strata=strata,shorten=.true.,info=info)
    call check(info==0.and.wt%has_truncation,'left truncation enabled')
    call check(all(wt%weight_trunc>=0.0_dp),'truncation weights nonnegative')
    do i=1,wt%n
        if(wt%count(i)==1)call check(abs(wt%weight_cens(i)-1.0_dp)<1e-12_dp,'normalized first censor weight')
    end do
    print '(a)','test_crprep: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
