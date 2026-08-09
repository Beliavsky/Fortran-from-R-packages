program test_example
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_options, lowrankqp_result, solve_low_rank_qp, &
        LRQP_CHOL, LRQP_LU, LRQP_SMW, LRQP_PFCF
    implicit none
    real(dp) :: vsq(6,6), vlo(6,3), d(6), a(3,6), b(3), u(6), target(6)
    type(lowrankqp_options) :: opt
    type(lowrankqp_result) :: res
    integer :: meth

    vsq=0.0_dp
    vsq(1,1)=1.0_dp; vsq(2,2)=1.0_dp; vsq(3,3)=1.0_dp
    vlo=0.0_dp
    vlo(1,1)=1.0_dp; vlo(2,2)=1.0_dp; vlo(3,3)=1.0_dp
    d=[0.0_dp,-5.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
    a=reshape([ -4.0_dp, 2.0_dp, 0.0_dp, &
                -3.0_dp, 1.0_dp,-2.0_dp, &
                 0.0_dp, 0.0_dp, 1.0_dp, &
                -1.0_dp, 0.0_dp, 0.0_dp, &
                 0.0_dp,-1.0_dp, 0.0_dp, &
                 0.0_dp, 0.0_dp,-1.0_dp ],shape(a))
    b=[-8.0_dp,2.0_dp,0.0_dp]
    u=100.0_dp
    target=[10.0_dp/21.0_dp,22.0_dp/21.0_dp,44.0_dp/21.0_dp,62.0_dp/21.0_dp,0.0_dp,0.0_dp]
    opt%tol=1.0e-9_dp
    opt%max_iter=200

    do meth=LRQP_LU,LRQP_CHOL
        opt%method=meth
        call solve_low_rank_qp(vsq,d,a,b,u,res,opt)
        call check(res,target,a,b)
    end do
    do meth=LRQP_SMW,LRQP_PFCF
        opt%method=meth
        call solve_low_rank_qp(vlo,d,a,b,u,res,opt)
        call check(res,target,a,b)
    end do
    print '(a)', 'PASS test_example'
contains
subroutine check(r,t,aa,bb)
    type(lowrankqp_result), intent(in) :: r
    real(dp), intent(in) :: t(:),aa(:,:),bb(:)
    if (.not.r%converged) error stop 'solver did not converge'
    if (maxval(abs(r%alpha-t)) > 2.0e-6_dp) error stop 'wrong alpha'
    if (maxval(abs(matmul(aa,r%alpha)-bb)) > 1.0e-7_dp) error stop 'equality residual'
    if (minval(r%alpha) < -1.0e-8_dp .or. maxval(r%alpha-100.0_dp) > 1.0e-8_dp) error stop 'bounds'
end subroutine check
end program test_example
