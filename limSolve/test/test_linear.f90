program test_linear
    use limsolve
    implicit none
    real(dp) :: dl(2),d(3),du(2),b(3,1),x(3,1),abd(3,3),afull(3,3)
    integer :: st
    dl=[-1.0_dp,-1.0_dp]; d=[2.0_dp,2.0_dp,2.0_dp]; du=dl
    b(:,1)=[1.0_dp,0.0_dp,1.0_dp]
    call solve_tridiag(dl,d,du,b,x,st)
    call check(st==LS_SUCCESS,'tridiag status')
    call check(maxval(abs(x(:,1)-[1.0_dp,1.0_dp,1.0_dp])) < 1.0e-12_dp,'tridiag')

    abd=0.0_dp
    abd(2,:)=[2.0_dp,2.0_dp,2.0_dp]
    abd(1,2:3)=-1.0_dp
    abd(3,1:2)=-1.0_dp
    call solve_banded(abd,1,1,b,x,st)
    call check(st==LS_SUCCESS,'banded status')
    call check(maxval(abs(x(:,1)-1.0_dp)) < 1.0e-12_dp,'banded')

    afull=0.0_dp
    afull(1,:)=[2.0_dp,-1.0_dp,0.0_dp]
    afull(2,:)=[-1.0_dp,2.0_dp,-1.0_dp]
    afull(3,:)=[0.0_dp,-1.0_dp,2.0_dp]
    call solve_banded(afull,1,1,b,x,st,full=.true.)
    call check(maxval(abs(x(:,1)-1.0_dp)) < 1.0e-12_dp,'banded full')
    print *, 'PASS test_linear'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_linear
