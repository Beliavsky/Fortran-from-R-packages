program test_relsurv_integration
    use mstate
    use relsurv_ratetable, only : ratetable_type, make_ratetable
    implicit none
    type(transition_map) :: tr,trn
    type(msdata_type) :: ms
    type(hazard_type) :: hz,hzn
    type(ratetable_type) :: tab
    integer,allocatable :: link(:,:),nlink(:)
    real(dp),allocatable :: xrate(:,:),pop(:,:),gt(:),hp(:),he(:),nr(:),ne(:),nc(:),se(:)
    integer :: info
    integer :: dims(1),fac(1),ncuts(1)
    real(dp) :: cuts(1,1),rates(1)

    call trans_comprisk(1,tr)
    ms%n=3
    allocate(ms%id(3),ms%from(3),ms%to(3),ms%trans(3),ms%status(3),ms%tstart(3),ms%tstop(3),ms%time(3))
    ms%id=[1,2,3];ms%from=1;ms%to=2;ms%trans=1;ms%status=1;ms%tstart=0.0_dp
    ms%tstop=[1.0_dp,2.0_dp,3.0_dp];ms%time=ms%tstop
    allocate(xrate(3,1));xrate=1.0_dp
    dims=1;fac=1;ncuts=0;cuts=0.0_dp;rates=0.01_dp
    tab=make_ratetable(dims,fac,cuts,ncuts,rates)
    call nelson_aalen_msdata(ms,tr,hz,info=info)
    call check(info==0,'Nelson-Aalen setup')
    call haz_function_relsurv(ms,xrate,1,tab,hz%time,gt,hp,he,nr,ne,nc,se,info)
    call check(info==0,'haz_function relsurv')
    call check(maxval(abs(hp-0.01_dp))<2e-12_dp,'constant population hazard increments')
    call msfit_relsurv(hz,ms,tr,[1],tab,xrate,hzn,trn,link,nlink,info,pop_haz_out=pop)
    call check(info==0,'msfit relsurv')
    call check(maxval(abs(hzn%haz(:,link(1,1))-[0.01_dp,0.02_dp,0.03_dp]))<2e-12_dp,'population split')
    call check(maxval(abs(hzn%haz(:,link(2,1))-(hz%haz(:,1)-[0.01_dp,0.02_dp,0.03_dp])))<2e-12_dp,'excess split')
    call check(maxval(abs(hzn%haz(:,link(1,1))+hzn%haz(:,link(2,1))-hz%haz(:,1)))<2e-12_dp,'split identity')
    print '(a)','test_relsurv_integration: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
