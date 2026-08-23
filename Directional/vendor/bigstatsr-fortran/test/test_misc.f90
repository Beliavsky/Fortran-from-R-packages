program test_misc
    use bigstatsr
    implicit none
    real(dp) :: b(2,3)
    real(dp), allocatable :: g(:)
    b=reshape([0.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,3.0_dp],[2,3])
    g=get_beta(b,'mean-wise')
    call check(maxval(abs(g-[1.0_dp,1.0_dp]))<1.0e-14_dp,'mean-wise')
    g=get_beta(b,'median-wise')
    call check(maxval(abs(g-[1.0_dp,0.0_dp]))<1.0e-14_dp,'median-wise')
    call check(block_size(1000,1,1.0_dp)>100000,'block_size')
    print *, 'test_misc: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok
        character(len=*),intent(in)::msg
        if(.not.ok) then
            print *, 'FAIL: ',trim(msg)
            error stop 1
        end if
    end subroutine check
end program test_misc
