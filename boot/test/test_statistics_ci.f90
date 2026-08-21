program test_statistics_ci
    use boot_kinds, only : dp
    use boot_statistics
    use boot_ci
    implicit none
    real(dp)::d(4,2),w(4),rho,l(4),lo(1),hi(1),conf(1),t(5),q(2),ranks(2)
    d(:,1)=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
    d(:,2)=2.0_dp*d(:,1)+1.0_dp
    w=1.0_dp
    rho=weighted_corr(d,w)
    if(abs(rho-1.0_dp)>1.0e-12_dp)error stop 1
    l=[-1.5_dp,-0.5_dp,0.5_dp,1.5_dp]
    if(abs(var_linear(l)-sum(l*l)/16.0_dp)>1.0e-12_dp)error stop 2
    if(abs(cum3(l,l,l,.false.))>1.0e-12_dp)error stop 3
    if(abs(inv_logit(logit(0.3_dp))-0.3_dp)>1.0e-12_dp)error stop 4
    t=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    call norm_inter(t,[0.5_dp,0.25_dp],q,ranks)
    if(abs(q(1)-3.0_dp)>1.0e-12_dp)error stop 5
    conf=[0.90_dp]
    call percentile_ci(t,conf,lo,hi)
    if(.not.(lo(1)<3.0_dp .and. hi(1)>3.0_dp))error stop 6
    print '(a)', 'test_statistics_ci: PASS'
end program test_statistics_ci
