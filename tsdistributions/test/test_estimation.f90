program test_estimation
  use tsdistributions
  implicit none
  type(rng_state)::rng
  type(distribution_parameters)::p
  type(parameter_specification)::spec
  type(distribution_fit)::fit
  real(dp),allocatable::x(:),v(:,:)
  integer::status
  call seed_rng(rng,918273_i8)
  p=distribution_parameters(mu=0.7_dp,sigma=1.4_dp)
  x=rdist('norm',600,rng,p)
  spec=distribution_modelspec(x,'norm')
  fit=estimate_distribution(x,spec,max_iterations=800)
  call assert_true(fit%status==tsd_success,'normal fit convergence')
  call assert_true(abs(fit%parameters%mu-p%mu)<0.15_dp,'normal mean recovery')
  call assert_true(abs(fit%parameters%sigma-p%sigma)<0.15_dp,'normal scale recovery')
  call assert_true(fit%aic<huge(1.0_dp).and.fit%bic<huge(1.0_dp),'information criteria')
  call information_covariance(fit,'QMLE',v,status)
  call assert_true(status==tsd_success.and.size(v,1)==2,'sandwich covariance')
  spec=distribution_modelspec(x,'std');spec%parameters%shape=8.0_dp;spec%estimate_shape=.false.
  fit=estimate_distribution(x,spec,max_iterations=600,use_hessian=.false.)
  call assert_true(fit%status==tsd_success,'fixed-shape Student fit')
  print '(a)','test_estimation: PASS'
contains
  subroutine assert_true(condition,message)
    logical,intent(in)::condition;character(len=*),intent(in)::message
    if(.not.condition)then;write(*,'(a)')'FAIL: '//message;error stop 1;end if
  end subroutine
end program
