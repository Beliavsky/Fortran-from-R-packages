program test_aft
  use survival
  implicit none
  real(dp)::time(5),x(5,1),logs(5),mu,sig
  integer::status(5)
  type(aft_result)::fit
  logs=[-1._dp,-0.5_dp,0._dp,0.5_dp,1._dp]
  time=exp(logs);x(:,1)=1._dp;status=1
  mu=sum(logs)/5._dp;sig=sqrt(sum((logs-mu)**2)/5._dp)
  call survreg_fit(time,status,x,'lognormal',fit,maxiter=25)
  if(abs(fit%coef(1)-mu)>2e-3_dp) error stop 'aft mean'
  if(abs(fit%scale-sig)>2e-3_dp) error stop 'aft scale'
  print *, 'test_aft PASS',fit%coef(1),fit%scale
end program
