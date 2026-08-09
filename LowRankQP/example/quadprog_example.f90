program quadprog_example
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_options, lowrankqp_result, solve_low_rank_qp, LRQP_CHOL
    implicit none
    real(dp) :: v(6,6), d(6), a(3,6), b(3), u(6)
    type(lowrankqp_options) :: opt
    type(lowrankqp_result) :: res
    v=0.0_dp; v(1,1)=1.0_dp; v(2,2)=1.0_dp; v(3,3)=1.0_dp
    d=[0.0_dp,-5.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
    a=reshape([-4.0_dp,2.0_dp,0.0_dp,-3.0_dp,1.0_dp,-2.0_dp,0.0_dp,0.0_dp,1.0_dp, &
        -1.0_dp,0.0_dp,0.0_dp,0.0_dp,-1.0_dp,0.0_dp,0.0_dp,0.0_dp,-1.0_dp],shape(a))
    b=[-8.0_dp,2.0_dp,0.0_dp]; u=100.0_dp
    opt%method=LRQP_CHOL; opt%verbose=.true.
    call solve_low_rank_qp(v,d,a,b,u,res,opt)
    print '(a,6(1x,f12.8))', 'alpha =',res%alpha
    print '(a,l1,a,i0)', 'converged = ',res%converged,' iterations = ',res%iterations
end program quadprog_example
