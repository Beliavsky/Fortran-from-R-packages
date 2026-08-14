program test_counting
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use survival
  implicit none
  real(dp) :: start(6),stop(6),x(6,1)
  integer :: status(6)
  type(coxph_result) :: fit
  type(survfit_result) :: sf
  start=[0._dp,0._dp,0._dp,1._dp,2._dp,3._dp]
  stop =[2._dp,3._dp,4._dp,4._dp,5._dp,6._dp]
  status=[1,0,1,0,1,0]
  x(:,1)=[1._dp,0._dp,1._dp,0._dp,1._dp,0._dp]
  call kaplan_meier_counting(start,stop,status,sf)
  if(size(sf%time)<3) error stop 'counting km'
  call coxph_fit_counting(start,stop,status,x,fit,'breslow',maxiter=20)
  if(.not.all(ieee_is_finite(fit%coef))) error stop 'counting cox finite'
  print *, 'test_counting PASS'
end program
