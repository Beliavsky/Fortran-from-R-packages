! SPDX-License-Identifier: GPL-3.0-only
program test_probabilities
  use svdnf
  use test_support
  implicit none
  real(dp) :: expected

  call assert_close(normal_pdf(0.0_dp),0.3989422804014327_dp,1.0e-14_dp,'normal density')
  call assert_close(normal_cdf(1.0_dp),0.8413447460685429_dp,1.0e-14_dp,'normal CDF')
  call assert_close(normal_quantile(0.975_dp),1.959963986120195_dp,2.0e-9_dp,'normal quantile')
  expected=1.0_dp-2.0_dp*exp(-1.0_dp)
  call assert_close(gamma_cdf(1.0_dp,2.0_dp,1.0_dp),expected,1.0e-13_dp,'gamma CDF')
  call assert_close(poisson_pmf(3,2.0_dp),exp(-2.0_dp)*8.0_dp/6.0_dp,1.0e-14_dp,'Poisson PMF')
  call assert_close(binomial_pmf(2,5,0.3_dp),0.3087_dp,1.0e-14_dp,'binomial PMF')

  write(*,'(a)') 'test_probabilities: PASS'
end program test_probabilities
