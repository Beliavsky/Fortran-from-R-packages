program test_relsurv_bootstrap
    use mstate
    use relsurv_ratetable, only : ratetable_type, make_ratetable
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    implicit none
    type(transition_map) :: tr,trn
    type(msdata_type) :: ms
    type(hazard_type) :: hz,hzn
    type(probtrans_type) :: pt,pt0
    type(relative_bootstrap_type) :: boot
    type(ratetable_type) :: tab
    integer,allocatable :: link(:,:),nlink(:),smap(:,:)
    real(dp),allocatable :: xrate(:,:)
    integer :: info,qpop,qexc
    integer :: dims(1),fac(1),ncuts(1)
    real(dp) :: cuts(1,1),rates(1),times(4)

    call trans_comprisk(1,tr)
    ms%n=8
    allocate(ms%id(8),ms%from(8),ms%to(8),ms%trans(8),ms%status(8),ms%tstart(8),ms%tstop(8),ms%time(8))
    ms%id=[1,2,3,4,5,6,7,8];ms%from=1;ms%to=2;ms%trans=1;ms%tstart=0.0_dp
    ms%status=[1,1,1,1,0,0,0,0]
    ms%tstop=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,2.5_dp,3.5_dp,4.5_dp,5.0_dp];ms%time=ms%tstop
    allocate(xrate(ms%n,1));xrate=1.0_dp
    dims=1;fac=1;ncuts=0;cuts=0.0_dp;rates=0.01_dp
    tab=make_ratetable(dims,fac,cuts,ncuts,rates)
    times=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    call modify_transition_relative(tr,[1],trn,link,nlink,smap,info)
    call check(info==0,'relative transition')
    qpop=link(1,1);qexc=link(2,1)
    call nelson_aalen_msdata(ms,tr,hz,info=info)
    call check(info==0,'point Nelson-Aalen')
    call msfit_relsurv(hz,ms,tr,[1],tab,xrate,hzn,trn,link,nlink,info)
    call check(info==0,'point relative msfit')
    call msboot_relsurv(ms,tr,[1],tab,xrate,times,24,boot,info,seed=7291,boot_original=.true.)
    call check(info==0,'relative bootstrap')
    call check(boot%nvalid(qpop)==24.and.boot%nvalid(qexc)==24,'all bootstrap samples retained')
    call check(all(ieee_is_finite(boot%varhaz(:,qpop))),'finite population bootstrap variance')
    call check(maxval(abs(boot%haz(:,qpop,:)-spread(0.01_dp*times,2,24)))<5e-12_dp, &
               'constant population hazards across bootstrap samples')
    call check(maxval(abs(boot%varhaz(:,qpop)))<1e-22_dp,'zero population bootstrap variance for constant rate')
    call check(maxval(boot%varhaz(:,qexc))>1e-6_dp,'positive excess bootstrap variance')
    call check(boot%has_original.and.boot%original_nvalid(1)==24,'original-hazard bootstrap retained')
    call check(maxval(boot%original_varhaz(:,1))>1e-6_dp,'original-hazard bootstrap variance')
    call check(maxval(abs(boot%haz(:,qpop,:)+boot%haz(:,qexc,:)-boot%original_haz(:,1,:)))<5e-12_dp, &
               'bootstrap split identity')
    call probtrans(hzn,trn,0.0_dp,pt0,variance=.false.,info=info)
    call check(info==0,'point probtrans')
    call probtrans_bootstrap(hzn,trn,boot,0.0_dp,pt,info=info)
    call check(info==0,'bootstrap probtrans')
    call check(maxval(abs(pt%p-pt0%p))<1e-12_dp,'bootstrap probtrans keeps point estimate')
    call check(maxval(pt%se(:,1,:))>1e-5_dp,'bootstrap transition-probability uncertainty')
    call check(maxval(abs(pt%se(:,2:3,:)))<1e-12_dp,'absorbing-state bootstrap standard errors')
    print '(a)','test_relsurv_bootstrap: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
