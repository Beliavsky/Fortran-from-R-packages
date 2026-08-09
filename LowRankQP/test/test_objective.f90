program test_objective
    use lowrankqp_kinds, only : dp
    use lowrankqp, only : lowrankqp_objective
    implicit none
    real(dp) :: v1(2,2), v2(3,2), d1(2), d2(3), x1(2), x2(3), f
    v1=reshape([2.0_dp,0.0_dp,0.0_dp,4.0_dp],shape(v1))
    d1=[-2.0_dp,-8.0_dp]; x1=[1.0_dp,2.0_dp]
    f=lowrankqp_objective(v1,d1,x1)
    if (abs(f+9.0_dp) > 1.0e-12_dp) error stop 'square objective'
    v2=0.0_dp; v2(1,1)=1.0_dp; v2(2,2)=2.0_dp
    d2=[1.0_dp,-1.0_dp,0.5_dp]; x2=[2.0_dp,1.0_dp,3.0_dp]
    f=lowrankqp_objective(v2,d2,x2)
    if (abs(f-6.5_dp) > 1.0e-12_dp) error stop 'low-rank objective'
    print '(a)', 'PASS test_objective'
end program test_objective
