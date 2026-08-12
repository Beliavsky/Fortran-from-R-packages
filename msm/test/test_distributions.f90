program test_distributions
    use msm, only : dp,dpexp,ppexp,qpexp,d2phase,p2phase,q2phase,h2phase,dtnorm,ptnorm,qtnorm
    implicit none
    real(dp) :: rate(3),tt(3),x,p,eps,fd,l1,mu1,mu2
    rate=[0.2_dp,0.5_dp,0.1_dp]; tt=[0.0_dp,2.0_dp,5.0_dp]; x=3.4_dp
    p=ppexp(x,rate,tt); call check(abs(qpexp(p,rate,tt)-x)<2e-12_dp,"piecewise exponential quantile")
    call check(abs(dpexp(1.0_dp,rate,tt)-0.2_dp*exp(-0.2_dp))<2e-14_dp,"piecewise exponential density")
    l1=0.3_dp; mu1=0.2_dp; mu2=0.7_dp; x=1.6_dp; eps=1e-6_dp
    fd=(p2phase(x+eps,l1,mu1,mu2)-p2phase(x-eps,l1,mu1,mu2))/(2.0_dp*eps)
    call check(abs(fd-d2phase(x,l1,mu1,mu2))<2e-10_dp,"two-phase density/CDF")
    p=p2phase(x,l1,mu1,mu2); call check(abs(q2phase(p,l1,mu1,mu2)-x)<2e-10_dp,"two-phase quantile")
    call check(abs(h2phase(x,l1,mu1,mu2)-d2phase(x,l1,mu1,mu2)/(1.0_dp-p))<2e-14_dp,"two-phase hazard")
    x=0.7_dp; p=ptnorm(x,0.0_dp,1.0_dp,-1.0_dp,2.0_dp)
    call check(abs(qtnorm(p,0.0_dp,1.0_dp,-1.0_dp,2.0_dp)-x)<2e-12_dp,"truncated normal quantile")
    fd=(ptnorm(x+eps,0.0_dp,1.0_dp,-1.0_dp,2.0_dp)-ptnorm(x-eps,0.0_dp,1.0_dp,-1.0_dp,2.0_dp))/(2.0_dp*eps)
    call check(abs(fd-dtnorm(x,0.0_dp,1.0_dp,-1.0_dp,2.0_dp))<2e-10_dp,"truncated normal density/CDF")
    print '(a)', 'test_distributions: PASS'
contains
    subroutine check(ok,msg)
        logical,intent(in)::ok; character(*),intent(in)::msg
        if(.not.ok) then; write(*,'(a)') 'FAIL: '//msg; error stop 1; end if
    end subroutine check
end program test_distributions
