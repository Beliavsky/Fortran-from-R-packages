program test_inference
  use dirichletreg, only : dp, dirichletreg_model, wald_confint, likelihood_ratio_test, &
       standardized_residuals, composite_residuals
  implicit none
  type(dirichletreg_model) :: m
  real(dp) :: lo(2),hi(2),dev,df,p,res(2,2),cres(2),y(2,2)
  integer :: stat,failures

  failures=0
  m%npar=2
  allocate(m%coefficients(2),m%se(2),m%mu(2,2),m%phi(2))
  m%coefficients=[1.0_dp,-0.5_dp]
  m%se=[0.2_dp,0.1_dp]
  m%mu=reshape([0.4_dp,0.6_dp,0.6_dp,0.4_dp],[2,2])
  m%phi=[9.0_dp,9.0_dp]
  y=reshape([0.45_dp,0.55_dp,0.55_dp,0.45_dp],[2,2])

  call wald_confint(m,0.95_dp,lo,hi,stat=stat)
  if(stat/=0 .or. abs(lo(1)-0.608007203091989_dp)>2.0e-12_dp .or. &
     abs(hi(1)-1.391992796908011_dp)>2.0e-12_dp) failures=failures+1

  call likelihood_ratio_test(-10.0_dp,5,-11.92_dp,4,dev,df,p,stat)
  if(stat/=0 .or. abs(dev-3.84_dp)>1.0e-12_dp .or. abs(df-1.0_dp)>1.0e-12_dp .or. &
     abs(p-0.0500435212487051_dp)>2.0e-10_dp) failures=failures+1

  call standardized_residuals(y,m,res,stat)
  call composite_residuals(y,m,cres,stat)
  if(stat/=0 .or. any(cres<0.0_dp)) failures=failures+1

  if(failures==0) then
    print '(a)', 'test_inference: PASS'
  else
    print '(a,i0)', 'test_inference: FAIL ',failures
    error stop 1
  end if
end program test_inference
