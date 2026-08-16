program test_distributions
  use fitdistrplus
  use test_support
  implicit none
  type(distribution_model) :: dist
  real(dp) :: par2(2),par1(1),x,p
  integer :: status

  call make_normal(dist);par2=[1.0_dp,2.0_dp]
  call assert_close(exp(dist%logpdf(1.0_dp,par2)),1.0_dp/(2.0_dp*sqrt(2.0_dp*pi_dp)),1.0e-12_dp,"normal density")
  call assert_close(dist%cdf(1.0_dp,par2),0.5_dp,1.0e-14_dp,"normal cdf")
  call assert_close(dist%quantile(0.5_dp,par2),1.0_dp,1.0e-14_dp,"normal quantile")
  call assert_close(dist%raw_moment(2,par2),5.0_dp,1.0e-14_dp,"normal moment")

  call make_gamma(dist);par2=[2.0_dp,3.0_dp];x=1.2_dp
  p=dist%cdf(x,par2)
  call assert_close(dist%cdf(dist%quantile(p,par2),par2),p,1.0e-9_dp,"gamma inverse cdf")
  call assert_close(dist%raw_moment(1,par2),2.0_dp/3.0_dp,1.0e-12_dp,"gamma mean")

  call make_weibull(dist);par2=[1.5_dp,2.0_dp];p=0.7_dp
  call assert_close(dist%cdf(dist%quantile(p,par2),par2),p,1.0e-12_dp,"weibull inverse cdf")

  call make_beta(dist);par2=[2.0_dp,5.0_dp];p=0.3_dp
  call assert_close(dist%cdf(dist%quantile(p,par2),par2),p,1.0e-9_dp,"beta inverse cdf")

  call make_poisson(dist);par1=[3.0_dp]
  call assert_close(dist%cdf(2.0_dp,par1),0.4231900811268435_dp,1.0e-10_dp,"poisson cdf")
  call assert_close(dist%raw_moment(2,par1),12.0_dp,1.0e-14_dp,"poisson moment")

  call make_negative_binomial(dist);par2=[4.0_dp,6.0_dp]
  call assert_close(dist%raw_moment(1,par2),6.0_dp,1.0e-14_dp,"negative binomial mean")
  call make_distribution("lnorm",dist,status)
  call assert_true(status==0 .and. trim(dist%name)=="lognormal","distribution factory")

  write(*,'(a)')"test_distributions: PASS"
end program test_distributions
