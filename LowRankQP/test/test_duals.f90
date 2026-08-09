program test_duals
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_options, lowrankqp_result, solve_low_rank_qp, LRQP_CHOL
    implicit none
    real(dp) :: v(2,2),d(2),a(1,2),b(1),u(2),station(2)
    type(lowrankqp_options) :: opt
    type(lowrankqp_result) :: res
    v=0.0_dp; v(1,1)=1.0_dp; v(2,2)=1.0_dp
    d=[-1.0_dp,-2.0_dp]; a(1,:)=[1.0_dp,1.0_dp]; b=1.0_dp; u=10.0_dp
    opt%method=LRQP_CHOL; opt%tol=1.0e-10_dp
    call solve_low_rank_qp(v,d,a,b,u,res,opt)
    if (.not.res%converged) error stop 'dual test convergence'
    station=matmul(v,res%alpha)+d+matmul(transpose(a),res%beta)+res%xi-res%zeta
    if (maxval(abs(station)) > 1.0e-7_dp) error stop 'stationarity'
    if (abs(sum(res%alpha)-1.0_dp) > 1.0e-8_dp) error stop 'equality'
    if (minval(res%xi) < -1.0e-10_dp .or. minval(res%zeta) < -1.0e-10_dp) error stop 'dual positivity'
    print '(a)', 'PASS test_duals'
end program test_duals
