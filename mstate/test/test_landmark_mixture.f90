program test_landmark_mixture
    use mstate
    implicit none
    type(transition_map)::tr
    type(msdata_type)::ms
    type(probtrans_type)::pt
    integer::mat(3,3),info
    mat=0;mat(1,3)=1;mat(2,3)=2
    call transition_from_matrix(mat,tr,info)
    ms%n=4
    allocate(ms%id(4),ms%from(4),ms%to(4),ms%trans(4),ms%status(4),ms%tstart(4),ms%tstop(4),ms%time(4))
    ms%id=[1,2,3,4];ms%from=[1,1,2,2];ms%to=3;ms%trans=[1,1,2,2];ms%status=1
    ms%tstart=0.0_dp;ms%tstop=[2.0_dp,3.0_dp,2.5_dp,4.0_dp];ms%time=ms%tstop
    call landmark_aj(ms,tr,1.0_dp,[1,2],pt,info)
    call check(info==0.and.pt%nt>1,'multi-start landmark fit')
    call check(maxval(abs(pt%p(1,1,:)-[0.5_dp,0.5_dp,0.0_dp]))<1e-12_dp,'empirical start mixture')
    call check(abs(pt%se(1,1,1)-0.25_dp)<1e-12_dp.and.abs(pt%se(1,1,2)-0.25_dp)<1e-12_dp, &
               'multinomial landmark SE')
    call check(all(pt%se(:,1,:)>=0.0_dp),'mixture SE nonnegative')
    print '(a)','test_landmark_mixture: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok;character(len=*),intent(in)::msg
        if(.not.ok)then;write(*,'(a,1x,a)')'FAIL:',msg;error stop 1;end if
    end subroutine
end program
