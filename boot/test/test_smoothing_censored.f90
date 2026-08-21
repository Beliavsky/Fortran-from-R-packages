program test_smoothing_censored
    use boot_kinds, only : dp
    use boot_smoothing
    use boot_censored
    implicit none
    integer::freq(4,3),strata(3),event(2)
    real(dp)::t(4),out(1,3),draws(20,2),ct(20,2)
    freq=reshape([1,1,1,1, 1,1,1,1, 1,1,1,1],[4,3])
    t=[0.0_dp,1.0_dp,2.0_dp,3.0_dp]
    strata=1
    call smooth_frequencies([1.5_dp],t,freq,strata,0.5_dp,out)
    if(abs(sum(out(1,:))-1.0_dp)>1.0e-12_dp)error stop 1
    call sample_product_limit([1.0_dp,2.0_dp],[0.5_dp,0.0_dp],20,2,draws)
    if(any(draws<1.0_dp).or.any(draws>2.0_dp))error stop 2
    event=[0,1]
    call sample_conditional_censoring([1.5_dp,0.5_dp],event,[1.0_dp,2.0_dp],[0.5_dp,0.0_dp],20,ct)
    if(any(abs(ct(:,1)-1.5_dp)>1.0e-12_dp))error stop 3
    if(any(ct(:,2)<1.0_dp).or.any(ct(:,2)>2.0_dp))error stop 4
    print '(a)', 'test_smoothing_censored: PASS'
end program test_smoothing_censored
