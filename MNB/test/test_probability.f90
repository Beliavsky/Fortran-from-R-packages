program test_probability
  use mnb, only : dp,nb_total_pmf,randomized_quantile_residuals,set_mnb_seed
  implicit none
  real(dp)::s,m,phi,mu,x(6,1),y(6),par(2),rq(2)
  integer::k
  phi=2.3_dp;mu=4.2_dp;s=0.0_dp;m=0.0_dp
  do k=0,100;s=s+nb_total_pmf(k,phi,mu);m=m+real(k,dp)*nb_total_pmf(k,phi,mu);end do
  if(abs(s-1.0_dp)>1.0e-9_dp)error stop 1
  if(abs(m-mu)>1.0e-8_dp)error stop 2
  x(:,1)=1.0_dp;y=[0.0_dp,1.0_dp,0.0_dp,2.0_dp,1.0_dp,0.0_dp];par=[2.0_dp,log(1.0_dp)]
  call set_mnb_seed(7);call randomized_quantile_residuals(par,y,x,2,3,rq)
  if(any(abs(rq)>10.0_dp))error stop 3
  print *, 'test_probability: PASS'
end program
