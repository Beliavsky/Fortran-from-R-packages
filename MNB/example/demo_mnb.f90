program demo_mnb
  use mnb, only : dp,set_mnb_seed,simulate_mnb,fit_mnb,mnb_fit_result
  implicit none
  integer,parameter::n=100,mi=3,nn=n*mi
  real(dp)::x(nn,2),truth(3),start(3)
  integer::y(nn),i
  type(mnb_fit_result)::fit
  x(:,1)=1.0_dp;do i=1,nn;x(i,2)=real(mod((i-1)/mi,2),dp);end do
  truth=[3.0_dp,0.3_dp,-0.4_dp];start=[2.0_dp,0.0_dp,0.0_dp]
  call set_mnb_seed(2026);call simulate_mnb(n,mi,x,truth,y);fit=fit_mnb(start,real(y,dp),x,n,mi)
  print '(a,3f11.5)', 'estimate:',fit%par
  print '(a,f14.6)', 'logLik:  ',fit%loglik
end program
