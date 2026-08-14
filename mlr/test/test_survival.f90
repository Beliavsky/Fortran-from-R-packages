program test_survival
  use mlr_kinds, only : dp
  use mlr_survival
  implicit none
  integer,parameter::n=60
  real(dp)::x(n,1),time(n)
  integer::status(n),i
  type(mlr_cox_model)::m
  real(dp),allocatable::risk(:)
  do i=1,n
    x(i,1)=real(mod(7*i,23),dp)/22.0_dp
    time(i)=exp(-0.8_dp*x(i,1))*(1.0_dp+0.20_dp*real(mod(11*i,17),dp)/16.0_dp)
    status(i)=merge(0,1,mod(i,7)==0)
  end do
  call fit_cox_learner(time,status,x,m)
  call predict_cox_risk(m,x,risk)
  if(.not.m%fit%converged)error stop 'cox convergence'
  if(measure_cindex(time,status,risk)<0.90_dp)error stop 'cindex'
  print *, 'test_survival: PASS'
end program test_survival
