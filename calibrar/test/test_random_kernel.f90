program test_random_kernel
  use calibrar, only : dp, set_seed, rtnorm_sample, dmvnorm_pdf, gaussian_kernel_2d, pi_dp
  implicit none
  real(dp) :: x,mu(2),sig(2,2),pt(2),v,xg(5),yg(5),z(5,5)
  integer :: i
  call set_seed(123)
  do i=1,100
    x=rtnorm_sample(0.0_dp,1.0_dp,-1.0_dp,2.0_dp)
    if(x < -1.0_dp .or. x > 2.0_dp) error stop "truncated normal bounds failed"
  end do
  mu=0.0_dp;sig=0.0_dp;sig(1,1)=1.0_dp;sig(2,2)=1.0_dp;pt=0.0_dp
  v=dmvnorm_pdf(pt,mu,sig)
  if(abs(v-1.0_dp/(2.0_dp*pi_dp))>1.0e-12_dp) error stop "dmvnorm failed"
  call gaussian_kernel_2d(mu,sig,[-1.0_dp,-1.0_dp],[1.0_dp,1.0_dp],5,5,xg,yg,z)
  if(abs(z(3,3)-v)>1.0e-12_dp) error stop "gaussian kernel failed"
  print *, "PASS test_random_kernel"
end program test_random_kernel
