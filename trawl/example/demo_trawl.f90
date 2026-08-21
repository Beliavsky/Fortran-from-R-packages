program demo_trawl
  use trawl
  implicit none
  type(trawl_spec)::spec
  type(trawl_fit_result)::fit
  integer,allocatable::x(:)
  real(dp),allocatable::xr(:)
  integer::status
  spec%kind='Exp';spec%lambda1=0.7_dp
  call set_trawl_seed(12345)
  call sim_univariate_trawl(1000.0_dp,x,burnin=20.0_dp,marginal='Poi',spec=spec,v=2.0_dp,status=status)
  if(status/=trawl_ok) error stop 'simulation failed'
  allocate(xr(size(x)));xr=real(x,dp)
  fit=fit_exptrawl(xr)
  print '(a,f10.5)','fitted exponential lambda = ',fit%lambda1
  print '(a,f10.5)','fitted trawl measure      = ',fit%lm
end program
