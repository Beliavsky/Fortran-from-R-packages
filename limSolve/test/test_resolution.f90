program test_resolution
    use limsolve
    implicit none
    type(resolution_result) :: r
    real(dp) :: a(3,2)
    a=0.0_dp; a(1,1)=1.0_dp; a(2,2)=1.0_dp
    r=resolution(a)
    call check(r%nsolvable==2,'rank')
    call check(maxval(abs(r%col-1.0_dp))<1.0e-10_dp,'column resolution')
    call check(maxval(abs(r%row-[1.0_dp,1.0_dp,0.0_dp]))<1.0e-10_dp,'row resolution')
    print *, 'PASS test_resolution'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_resolution
