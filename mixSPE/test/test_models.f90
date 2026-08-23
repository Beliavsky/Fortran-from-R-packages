program test_models
  use mixspe, only : dp, rpe, spe_model, em_fit, emgr_fit
  implicit none
  character(len=5), parameter :: models(8) = [character(len=5) :: &
    'EIIE ', 'VIIE ', 'EEIE ', 'VVIE ', 'EEEE ', 'EEVE ', 'VVEE ', 'VVVE ']
  character(len=5), parameter :: skew_models(2) = [character(len=5) :: 'EIIVV','EEEVV']
  real(dp), allocatable :: a(:,:),b(:,:),x(:,:),bics(:,:)
  real(dp) :: s(2,2)
  type(spe_model) :: fit,best
  integer :: i
  s=0.0_dp; s(1,1)=1.0_dp; s(2,2)=0.7_dp
  a=rpe(80,1.2_dp,[0.0_dp,0.0_dp],s)
  b=rpe(80,0.8_dp,[2.8_dp,0.3_dp],s)
  allocate(x(160,2)); x(1:80,:)=a; x(81:160,:)=b
  do i=1,size(models)
    call em_fit(x,2,trim(models(i)),fit,max_iter=35,tol=1.0e-4_dp)
    call check(fit%loglik>-huge(1.0_dp)/100.0_dp,trim(models(i)))
    call check(all(fit%lam>0.0_dp),'positive eigenvalues '//trim(models(i)))
  end do
  do i=1,size(skew_models)
    call em_fit(x,2,trim(skew_models(i)),fit,max_iter=20,tol=1.0e-4_dp)
    call check(fit%loglik>-huge(1.0_dp)/100.0_dp,trim(skew_models(i)))
  end do
  call emgr_fit(x,[1,2],models,best,bics,max_iter=25,tol=1.0e-4_dp)
  call check(size(bics,1)==8 .and. size(bics,2)==2,'emgr dimensions')
  call check(best%g>=1 .and. best%g<=2,'emgr best model')
  print '(a)', 'test_models: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then; print '(a)', 'FAIL: '//trim(msg); error stop 1; end if
  end subroutine
end program
