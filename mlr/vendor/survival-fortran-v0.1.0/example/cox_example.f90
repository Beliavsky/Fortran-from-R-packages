program cox_example
  use survival
  implicit none
  real(dp) :: time(8),x(8,1)
  integer :: status(8)
  type(coxph_result) :: fit
  time=[1._dp,2._dp,3._dp,4._dp,5._dp,6._dp,7._dp,8._dp]
  status=[1,1,1,1,1,1,0,0]
  x(:,1)=[1._dp,1._dp,1._dp,1._dp,0._dp,0._dp,0._dp,0._dp]
  call coxph_fit(time,status,x,fit,'efron',maxiter=60)
  print '(a,l1)', 'converged: ',fit%converged
  print '(a,f12.6)', 'coefficient: ',fit%coef(1)
  print '(a,f12.6)', 'log likelihood: ',fit%loglik
end program
