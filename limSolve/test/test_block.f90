program test_block
    use limsolve
    implicit none
    real(dp) :: top(1,2),blocks(1,3,2),bot(1,2),b(4,1),x(4,1)
    integer :: st
    ! Full matrix assembled by solve_block is tridiagonal 4x4 with diagonal 2 and off-diagonal -1.
    top(1,:)=[2.0_dp,-1.0_dp]
    blocks(:,:,1)=reshape([-1.0_dp,2.0_dp,-1.0_dp],[1,3])
    blocks(:,:,2)=reshape([-1.0_dp,2.0_dp,-1.0_dp],[1,3])
    bot(1,:)=[-1.0_dp,2.0_dp]
    b(:,1)=[1.0_dp,0.0_dp,0.0_dp,1.0_dp]
    call solve_block(top,blocks,bot,b,2,x,st)
    call check(st==LS_SUCCESS,'block status')
    call check(maxval(abs(x(:,1)-1.0_dp))<1.0e-12_dp,'block solution')
    print *, 'PASS test_block'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_block
