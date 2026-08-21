program test_saddle
    use boot_kinds, only : dp, pi
    use boot_saddle
    implicit none
    real(dp)::a(4,1),u(1),mu(4),pdf,cdf,zeta(1),expected
    integer::info
    a(:,1)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    u=[6.0_dp]
    mu=1.0_dp
    call multinomial_saddlepoint(a,u,mu,pdf,cdf,zeta,info=info)
    if(info/=0)error stop 1
    expected=1.0_dp/sqrt(2.0_dp*pi*5.0_dp)
    if(abs(pdf-expected)>1.0e-8_dp)error stop 2
    if(abs(cdf-0.5_dp)>1.0e-8_dp)error stop 3
    if(abs(zeta(1))>1.0e-8_dp)error stop 4
    print '(a)', 'test_saddle: PASS'
end program test_saddle
