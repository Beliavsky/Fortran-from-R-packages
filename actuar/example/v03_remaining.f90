module v03_example_callbacks
    use actuar_kinds, only: dp
    implicit none
contains
    function exp_cdf(x,par) result(v)
        real(dp),intent(in)::x,par(:)
        real(dp)::v
        if(x<=0.0_dp)then;v=0.0_dp;else;v=1.0_dp-exp(-par(1)*x);end if
    end function exp_cdf
    function exp_pdf(x,par) result(v)
        real(dp),intent(in)::x,par(:)
        real(dp)::v
        if(x<0.0_dp)then;v=0.0_dp;else;v=par(1)*exp(-par(1)*x);end if
    end function exp_pdf
end module v03_example_callbacks

program v03_remaining
    use actuar
    use v03_example_callbacks
    implicit none
    type(mde_result_t)::fit
    type(coverage_spec_t)::cov
    real(dp)::bounds(5),counts(4),par(1),payment
    bounds=[0.0_dp,-log(0.75_dp)/2.0_dp,-log(0.5_dp)/2.0_dp,-log(0.25_dp)/2.0_dp,10.0_dp]
    counts=25.0_dp
    fit=mde_grouped_chisq(bounds,counts,exp_cdf,[1.0_dp],lower=[0.1_dp],upper=[5.0_dp])
    print '(a,f10.6)','MDE exponential rate: ',fit%estimate(1)
    cov%deductible=1.0_dp;cov%limit=3.0_dp;cov%coinsurance=0.8_dp;cov%inflation=0.1_dp;cov%per_loss=.true.
    par=[2.0_dp];payment=0.6_dp
    print '(a,f10.6)','Covered-loss CDF at 0.6: ',coverage_cdf(payment,par,exp_cdf,cov)
    print '(a,f10.6)','Covered-loss density at 0.6: ',coverage_pdf(payment,par,exp_pdf,exp_cdf,cov)
    print '(a,f10.6)','Probability mass at payment limit: ', &
        coverage_pdf(cov%maximum_payment(),par,exp_pdf,exp_cdf,cov)
end program v03_remaining
