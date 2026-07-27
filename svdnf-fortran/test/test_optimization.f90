! SPDX-License-Identifier: GPL-3.0-only
program test_optimization
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use svdnf
  use test_support
  implicit none
  type(svm_dynamics) :: truth, starting
  type(simulation_result) :: simulated
  type(filter_result) :: initial_filter
  type(optimization_result) :: fit
  real(dp) :: initial_loglik

  truth=dynamics_svm('Taylor',phi=0.90_dp,theta=-1.0_dp,sigma=0.16_dp)
  simulated=model_simulate(truth,45,initial_volatility=-1.0_dp,seed=4815)
  starting=dynamics_svm('Taylor',phi=0.75_dp,theta=-0.7_dp,sigma=0.28_dp)
  initial_filter=dnf_filter(starting,simulated%returns,n=10)
  call assert_true(initial_filter%ok,'initial optimization filter')
  initial_loglik=initial_filter%log_likelihood
  fit=dnf_optimize(starting,simulated%returns,n=10,max_iterations=45, &
    tolerance=2.0e-5_dp,calculate_hessian=.true.)
  call assert_true(fit%ok,'optimization result')
  call assert_true(ieee_is_finite(fit%log_likelihood),'finite fitted likelihood')
  call assert_true(fit%log_likelihood>=initial_loglik-1.0e-8_dp,'optimizer improves likelihood')
  call assert_true(size(fit%parameters)==3,'Taylor parameter count')
  call assert_true(fit%parameters(1)>-1.0_dp .and. fit%parameters(1)<1.0_dp,'valid fitted phi')
  call assert_true(fit%parameters(3)>0.0_dp,'valid fitted sigma')
  call assert_true(size(fit%hessian,1)==3,'Hessian size')
  call assert_close(maxval(abs(fit%hessian-transpose(fit%hessian))),0.0_dp,1.0e-12_dp,'Hessian symmetry')

  write(*,'(a)') 'test_optimization: PASS'
end program test_optimization
