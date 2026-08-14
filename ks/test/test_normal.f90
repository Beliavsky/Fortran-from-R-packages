! SPDX-License-Identifier: GPL-2.0-only
program test_normal
  use ks, only: dp, pi, normal_pdf, normal_cdf, normal_quantile, mvn_pdf, mvn_derivative_tensor, &
                dnorm_mixture, normal_mixture_moments
  implicit none
  real(dp)::p,q,f,mu2(2),S(2,2),x2(2),expect,mean(2),cov(2,2)
  real(dp),allocatable::d1(:),d2(:)
  real(dp)::mus1(2),sd1(2),pr(2),mums(2,2),sigs(2,2,2)
  integer::info
  p=0.123456789_dp;q=normal_quantile(p)
  if(abs(normal_cdf(q)-p)>2e-10_dp) error stop 'normal quantile/cdf'
  if(abs(normal_pdf(0.0_dp)-1.0_dp/sqrt(2.0_dp*pi))>1e-14_dp) error stop 'normal pdf'
  mu2=[0.2_dp,-0.4_dp];S=reshape([1.2_dp,0.3_dp,0.3_dp,0.8_dp],[2,2]);x2=[0.1_dp,0.7_dp]
  f=mvn_pdf(x2,mu2,S);expect=0.07098707029856734_dp
  if(abs(f-expect)>2e-14_dp) error stop 'mvn pdf'
  call mvn_derivative_tensor(x2,mu2,S,1,d1,info);if(info/=0) error stop 'mvn d1 info'
  call mvn_derivative_tensor(x2,mu2,S,2,d2,info);if(info/=0) error stop 'mvn d2 info'
  ! Finite difference first derivative in coordinate 1.
  if(abs(d1(1)-(mvn_pdf(x2+[1e-6_dp,0.0_dp],mu2,S)-mvn_pdf(x2-[1e-6_dp,0.0_dp],mu2,S))/2e-6_dp)>2e-8_dp) &
      error stop 'mvn derivative'
  mus1=[-1.0_dp,2.0_dp];sd1=[0.5_dp,1.2_dp];pr=[0.3_dp,0.7_dp]
  if(abs(dnorm_mixture(0.4_dp,mus1,sd1,pr) - &
     (pr(1)*normal_pdf(0.4_dp,mus1(1),sd1(1))+pr(2)*normal_pdf(0.4_dp,mus1(2),sd1(2))))>1e-14_dp) error stop 'mixture pdf'
  mums=reshape([-1.0_dp,0.2_dp,2.0_dp,1.5_dp],[2,2],order=[2,1])
  sigs=0.0_dp;sigs(:,:,1)=reshape([1.0_dp,0.2_dp,0.2_dp,0.7_dp],[2,2]);sigs(:,:,2)=reshape([0.5_dp,-0.1_dp,-0.1_dp,1.1_dp],[2,2])
  call normal_mixture_moments(mums,sigs,pr,mean,cov)
  if(maxval(abs(mean-(pr(1)*mums(1,:)+pr(2)*mums(2,:))))>1e-14_dp) error stop 'mixture moments mean'
  if(any([(cov(1,1)<=0.0_dp),(cov(2,2)<=0.0_dp)])) error stop 'mixture covariance'
  print *, 'test_normal: PASS'
end program
