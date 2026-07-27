! SPDX-License-Identifier: GPL-3.0-only
program estimate_taylor
  use svdnf
  implicit none
  type(svm_dynamics) :: truth, starting
  type(simulation_result) :: simulated
  type(optimization_result) :: fit
  character(len=24), allocatable :: names(:)
  integer :: i

  truth=dynamics_svm('Taylor',phi=0.94_dp,theta=-1.1_dp,sigma=0.17_dp)
  simulated=model_simulate(truth,100,initial_volatility=-1.1_dp,seed=2026)
  starting=dynamics_svm('Taylor',phi=0.8_dp,theta=-0.8_dp,sigma=0.25_dp)
  fit=dnf_optim(starting,simulated%returns,n=16,max_iterations=100,tolerance=1.0e-5_dp)
  if (.not. fit%ok) error stop trim(fit%message)
  names=model_parameter_names(starting)
  write(*,'(a,f14.6)') 'maximized log likelihood: ',fit%log_likelihood
  do i=1,size(fit%parameters)
    write(*,'(a24,2x,f14.7)') trim(names(i)),fit%parameters(i)
  end do
end program estimate_taylor
