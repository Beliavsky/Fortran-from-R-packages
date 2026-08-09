program test_inverse
    use limsolve
    implicit none
    type(solve_result) :: r
    real(dp) :: a(2,2), bmat(2,1), xmat(2,1)
    real(dp) :: ann(3,2), bnn(3)
    real(dp) :: g(2,2), h(2)

    a = reshape([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[2,2])
    bmat(:,1) = [1.0_dp,2.0_dp]
    call solve_generalized(a,bmat,xmat)
    call check(maxval(abs(xmat(:,1)-[0.2_dp,0.4_dp])) < 1.0e-7_dp,'generalized inverse')

    ann = reshape([1.0_dp,0.0_dp,1.0_dp, 0.0_dp,1.0_dp,1.0_dp],[3,2])
    bnn = [1.0_dp,-1.0_dp,0.0_dp]
    call nnls(ann,bnn,r)
    call check(r%succeeded(),'nnls status')
    call check(all(r%x >= -1.0e-10_dp),'nnls nonnegative')
    call check(maxval(abs(r%x-[0.5_dp,0.0_dp])) < 1.0e-6_dp,'nnls solution')

    g = 0.0_dp; g(1,1)=1.0_dp; g(2,2)=1.0_dp
    h = [1.0_dp,2.0_dp]
    call ldp(g,h,r)
    call check(r%succeeded(),'ldp status')
    call check(maxval(abs(r%x-[1.0_dp,2.0_dp])) < 1.0e-7_dp,'ldp solution')

    print *, 'PASS test_inverse'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(len=*),intent(in)::msg
        if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
    end subroutine check
end program test_inverse
