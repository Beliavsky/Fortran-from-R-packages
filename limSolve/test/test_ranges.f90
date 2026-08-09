program test_ranges
    use limsolve
    implicit none
    type(solve_result) :: r
    type(range_result) :: rr
    real(dp) :: e(1,2),f(1),g(2,2),h(2),c(2),eqa(1,2),vv(2,1)
    real(dp) :: samples(2,2)
    integer :: st
    e(1,:)=[1.0_dp,1.0_dp]; f=[1.0_dp]
    g=0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp; h=0.0_dp
    c=[1.0_dp,2.0_dp]
    call linp(e,f,g,h,c,r,ispos=.true.)
    call check(r%succeeded(),'linp status')
    call check(maxval(abs(r%x-[1.0_dp,0.0_dp])) < 1.0e-8_dp,'linp solution')

    rr=xranges(e,f,g,h,ispos=.true.,central=.true.)
    call check(rr%status==LS_SUCCESS,'xranges status')
    call check(maxval(abs(rr%range(:,1))) < 1.0e-8_dp,'xranges min')
    call check(maxval(abs(rr%range(:,2)-1.0_dp)) < 1.0e-8_dp,'xranges max')

    eqa(1,:)=[1.0_dp,-1.0_dp]
    rr=varranges(e,f,g,h,eqa,ispos=.true.)
    call check(abs(rr%range(1,1)+1.0_dp)<1.0e-8_dp,'varrange min')
    call check(abs(rr%range(1,2)-1.0_dp)<1.0e-8_dp,'varrange max')

    samples(1,:)=[1.0_dp,0.0_dp]; samples(2,:)=[0.25_dp,0.75_dp]
    call varsample(samples,eqa,vv,status=st)
    call check(st==LS_SUCCESS,'varsample status')
    call check(maxval(abs(vv(:,1)-[1.0_dp,-0.5_dp]))<1.0e-12_dp,'varsample')
    print *, 'PASS test_ranges'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_ranges
