program test_fit_sim
  use mnb, only : dp,mnb_loglik,fit_mnb,simulate_mnb,set_mnb_seed,mnb_fit_result
  implicit none
  integer,parameter::n=40,mi=3,p=2,nn=n*mi
  real(dp)::x(nn,p),par(3),start(3),ll
  integer::y(nn),i,j,k
  type(mnb_fit_result)::fit
  x(:,1)=1.0_dp
  do i=1,n
    do j=1,mi;k=(i-1)*mi+j;x(k,2)=merge(1.0_dp,0.0_dp,mod(i,2)==0);end do
  end do
  par=[2.5_dp,0.2_dp,0.35_dp];call set_mnb_seed(12345);call simulate_mnb(n,mi,x,par,y)
  ll=mnb_loglik(par,real(y,dp),x,n,mi);if(.not.(ll<0.0_dp))error stop 1
  start=[2.0_dp,0.0_dp,0.0_dp];fit=fit_mnb(start,real(y,dp),x,n,mi,maxit=300)
  if(fit%par(1)<=0.0_dp)error stop 2
  if(fit%loglik < mnb_loglik(start,real(y,dp),x,n,mi)-1.0e-6_dp)error stop 3
  if(any(fit%se<0.0_dp))error stop 4
  print *, 'test_fit_sim: PASS'
end program
