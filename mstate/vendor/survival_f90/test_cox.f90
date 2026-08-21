program test_cox
  use survival
  implicit none
  real(dp) :: time(8),x(8,1)
  integer :: status(8)
  type(coxph_result) :: fit
  time=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp]
  status=[1,1,1,1,1,1,0,0]
  x(:,1)=[1._dp,1._dp,1._dp,1._dp,0._dp,0._dp,0._dp,0._dp]
  call coxph_fit(time,status,x,fit,'efron',maxiter=60)
  if(.not.fit%converged) error stop 'cox convergence'
  if(fit%coef(1)<=0.0_dp) error stop 'cox sign'
  if(fit%loglik<fit%loglik_initial) error stop 'cox likelihood'
  print *, 'test_cox PASS', fit%coef(1)
end program
