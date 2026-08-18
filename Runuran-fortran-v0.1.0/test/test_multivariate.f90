! SPDX-License-Identifier: GPL-2.0-or-later
program test_multivariate
  use runuran
  implicit none
  type(rng_state)::rng
  type(multivariate_distribution)::d
  real(dp)::mu(2),cov(2,2),x(2),lp
  real(dp),allocatable::z(:,:)
  integer::i
  mu=[1.0_dp,-2.0_dp];cov=reshape([1.0_dp,0.4_dp,0.4_dp,2.0_dp],[2,2])
  d=udmultinormal(mu,cov);lp=d%logpdf(mu)
  if(abs(lp-(-log(2.0_dp*pi)-0.5_dp*log(1.84_dp)))>2e-12_dp)then;print *,'mv logpdf',lp;error stop 1;end if
  call rng_seed(rng,99_i8);allocate(z(2,6000));do i=1,size(z,2);call d%sample(rng,z(:,i));end do
  x=sum(z,dim=2)/real(size(z,2),dp)
  if(any(abs(x-mu)>0.09_dp))then;print *,'mv mean',x;error stop 1;end if
  print *,'test_multivariate: PASS'
end program
