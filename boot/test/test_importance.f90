program test_importance
    use boot_kinds, only : dp
    use boot_importance
    implicit none
    real(dp)::t(5),w(5),raw(2),rat(2),reg(2),pr(1),pp(1),pg(1),q1(1),q2(1),q3(1)
    t=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
    w=1.0_dp
    call importance_moments(t,w,raw,rat,reg)
    if(maxval(abs(raw-[3.0_dp,2.5_dp]))>1.0e-12_dp)error stop 1
    if(maxval(abs(rat-[3.0_dp,2.5_dp]))>1.0e-12_dp)error stop 2
    call importance_probability(t,w,[3.0_dp],pr,pp,pg)
    if(abs(pr(1)-0.6_dp)>1.0e-12_dp .or. abs(pp(1)-0.6_dp)>1.0e-12_dp)error stop 3
    call importance_quantile(t,w,[0.5_dp],q1,q2,q3)
    if(q2(1)<2.0_dp .or. q2(1)>3.0_dp)error stop 4
    print '(a)', 'test_importance: PASS'
end program test_importance
