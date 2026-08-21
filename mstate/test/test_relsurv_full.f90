program test_relsurv_full
    use mstate
    use relsurv_ratetable, only : ratetable_type, make_ratetable
    implicit none
    type(transition_map) :: tr
    type(msdata_type) :: ms
    type(hazard_type) :: hz,hza
    type(relative_msfit_type) :: fit
    type(ratetable_type) :: tab
    real(dp),allocatable :: xrate(:,:)
    integer :: info
    integer :: dims(1),fac(1),ncuts(1)
    real(dp) :: cuts(1,1),rates(1),expected(3)

    call trans_comprisk(1,tr)
    ms%n=3
    allocate(ms%id(3),ms%from(3),ms%to(3),ms%trans(3),ms%status(3),ms%tstart(3),ms%tstop(3),ms%time(3))
    ms%id=[1,2,3];ms%from=1;ms%to=2;ms%trans=1;ms%status=1;ms%tstart=0.0_dp
    ms%tstop=[1.0_dp,2.0_dp,3.0_dp];ms%time=ms%tstop
    allocate(xrate(ms%n,1));xrate=1.0_dp
    dims=1;fac=1;ncuts=0;cuts=0.0_dp;rates=1.0e-4_dp
    tab=make_ratetable(dims,fac,cuts,ncuts,rates)
    call nelson_aalen_msdata(ms,tr,hz,info=info)
    call check(info==0,'Nelson-Aalen setup')
    call hazard_add_times(hz,[1.5_dp],hza)
    call check(hza%nt==4.and.abs(hza%haz(2,1)-hz%haz(1,1))<1e-12_dp,'add-times carry-forward')
    call msfit_relsurv_full(hz,ms,tr,[1],tab,xrate,fit,variance_mode='fixed',time_format='years',info=info, &
                            add_times=[1.5_dp])
    call check(info==0,'full relative msfit')
    call check(fit%fit%nt==4,'full relative add-times grid')
    call check(abs(fit%population_haz(2,1)-1.0e-4_dp*365.241_dp*1.5_dp)<2e-11_dp, &
               'add-time population hazard')
    expected=1.0e-4_dp*365.241_dp*[1.0_dp,2.0_dp,3.0_dp]
    call check(maxval(abs(fit%population_haz([1,3,4],1)-expected))<2e-11_dp,'years-to-days population hazard scaling')
    call check(.not.fit%has_bootstrap.and.fit%variance_mode=='fixed','fixed variance mode')
    call msfit_relsurv_full(hz,ms,tr,[1],tab,xrate,fit,variance_mode='both',b=12,seed=31, &
                            time_format='years',info=info)
    call check(info==0.and.fit%has_bootstrap,'both variance mode')
    call check(all(fit%bootstrap%nvalid>=2),'bootstrap sample availability')
    call check(maxval(abs(fit%bootstrap_fit%haz-fit%fit%haz))<1e-12_dp,'bootstrap fit keeps point hazards')
    print '(a)','test_relsurv_full: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
