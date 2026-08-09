program test_unconstrained
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_options, lowrankqp_result, solve_low_rank_qp, &
        LRQP_LU, LRQP_CHOL, LRQP_SMW, LRQP_PFCF
    implicit none
    real(dp) :: vs(2,2), vl(3,2), ds(2), dl(3), us(2), ul(3)
    real(dp) :: a0s(0,2), a0l(0,3), b0(0)
    type(lowrankqp_options) :: opt
    type(lowrankqp_result) :: res
    integer :: meth

    vs=0.0_dp; vs(1,1)=2.0_dp; vs(2,2)=4.0_dp
    ds=[-2.0_dp,-8.0_dp]; us=10.0_dp
    opt%tol=1.0e-10_dp
    do meth=LRQP_LU,LRQP_CHOL
        opt%method=meth
        call solve_low_rank_qp(vs,ds,a0s,b0,us,res,opt)
        if (.not.res%converged) error stop 'square unconstrained did not converge'
        if (maxval(abs(res%alpha-[1.0_dp,2.0_dp])) > 2.0e-7_dp) error stop 'square unconstrained solution'
    end do

    vl=0.0_dp; vl(1,1)=sqrt(2.0_dp); vl(2,2)=2.0_dp
    dl=[-2.0_dp,-8.0_dp,1.0_dp]; ul=10.0_dp
    do meth=LRQP_SMW,LRQP_PFCF
        opt%method=meth
        call solve_low_rank_qp(vl,dl,a0l,b0,ul,res,opt)
        if (.not.res%converged) error stop 'low-rank unconstrained did not converge'
        if (maxval(abs(res%alpha-[1.0_dp,2.0_dp,0.0_dp])) > 2.0e-7_dp) error stop 'low-rank unconstrained solution'
    end do
    print '(a)', 'PASS test_unconstrained'
end program test_unconstrained
