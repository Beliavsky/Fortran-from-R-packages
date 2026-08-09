program test_status
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_options, lowrankqp_result, solve_low_rank_qp, LRQP_CHOL
    implicit none
    real(dp) :: v(2,2),d(2),a0(0,2),b0(0),ubad(1),u(2)
    type(lowrankqp_options) :: opt
    type(lowrankqp_result) :: res
    v=0.0_dp; v(1,1)=1.0_dp; v(2,2)=1.0_dp; d=0.0_dp; ubad=1.0_dp; u=1.0_dp
    call solve_low_rank_qp(v,d,a0,b0,ubad,res)
    if (res%status /= -1) error stop 'dimension status'
    opt%method=99
    call solve_low_rank_qp(v,d,a0,b0,u,res,opt)
    if (res%status /= -3) error stop 'method status'
    opt%method=LRQP_CHOL; opt%max_iter=1; opt%tol=1.0e-30_dp
    call solve_low_rank_qp(v,[-0.2_dp,-0.3_dp],a0,b0,u,res,opt)
    if (res%status /= 1 .or. res%converged) error stop 'iteration limit status'
    print '(a)', 'PASS test_status'
end program test_status
