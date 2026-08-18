program test_api
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  real(dp) :: alpha(3),d(100,3),mu(2),root(2,2),tt(50,2),v(2,2)
  type(wishart_result) :: wr
  call rng_seed(17)
  alpha=[1.0_dp,2.0_dp,3.0_dp]; d=rdirichlet(100,alpha)
  if (maxval(abs(sum(d,dim=2)-1.0_dp))>1.0e-12_dp) error stop "test_api: dirichlet"
  mu=0.0_dp; root=0.0_dp; root(1,1)=1.0_dp; root(2,2)=1.0_dp; tt=rmvst(50,7.0_dp,mu,root)
  if (.not.all(ieee_is_finite(tt))) error stop "test_api: rmvst"
  v=0.0_dp; v(1,1)=1.0_dp; v(2,2)=1.0_dp; wr=rwishart(6.0_dp,v)
  if (.not.all(ieee_is_finite(wr%w)) .or. wr%w(1,1)<=0.0_dp) error stop "test_api: wishart"
  print *, "test_api: PASS"
end program test_api
