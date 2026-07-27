! SPDX-License-Identifier: GPL-3.0-only
program test_models
  use svdnf
  use test_support
  use test_custom_callbacks
  implicit none
  type(svm_dynamics) :: dynamics, copied, custom
  type(simulation_result) :: custom_simulated
  type(filter_result) :: custom_filtered, custom_initial
  type(optimization_result) :: custom_fit
  type(grid_type) :: grids
  real(dp), allocatable :: parameters(:)
  logical :: ok
  character(len=160) :: message

  dynamics=dynamics_svm('Heston',mu=0.05_dp,kappa=2.0_dp,theta=0.04_dp, &
    sigma=0.3_dp,rho=-0.6_dp,h=1.0_dp/252.0_dp)
  call validate_dynamics(dynamics,ok,message)
  call assert_true(ok,'Heston validation')
  call assert_close(evaluate_mu_y(dynamics,0.04_dp), &
    (0.05_dp-0.02_dp)/252.0_dp,1.0e-14_dp,'Heston return drift')
  call assert_close(evaluate_sigma_y(dynamics,0.04_dp),sqrt(0.04_dp/252.0_dp), &
    1.0e-14_dp,'Heston return diffusion')
  call assert_close(evaluate_mu_x(dynamics,0.03_dp),0.03_dp+2.0_dp*(0.01_dp)/252.0_dp, &
    1.0e-14_dp,'Heston state drift')

  parameters=parameter_vector(dynamics)
  copied=dynamics_svm('Heston')
  call set_parameter_vector(copied,parameters,ok)
  call assert_true(ok,'parameter unpacking')
  call assert_vector_close(parameter_vector(copied),parameters,1.0e-14_dp,'parameter round trip')

  grids=grid_maker(dynamics,n=17,k=5,r=2)
  call assert_true(size(grids%var_mid_points)==17,'Heston grid size')
  call assert_true(size(grids%jump_counts)==1,'Heston no-jump count')
  call assert_true(all(grids%var_mid_points(2:)>grids%var_mid_points(:16)),'increasing grid')

  custom=dynamics_svm('Custom',rho=0.25_dp)
  call set_custom_dynamics(custom,custom_mu_y,custom_sigma_y,custom_mu_x,custom_sigma_x, &
    mu_y_parameters=[0.1_dp],sigma_y_parameters=[0.2_dp], &
    mu_x_parameters=[0.9_dp],sigma_x_parameters=[0.3_dp])
  call validate_dynamics(custom,ok,message)
  call assert_true(ok,'custom callback validation')
  call assert_close(evaluate_mu_y(custom,2.0_dp),0.3_dp,1.0e-14_dp,'custom mu_y')
  call assert_close(evaluate_mu_x(custom,2.0_dp),1.8_dp,1.0e-14_dp,'custom mu_x')

  custom_simulated=model_simulate(custom,25,initial_volatility=0.0_dp,seed=31415)
  call assert_true(custom_simulated%ok,'custom model simulation')
  custom_filtered=dnf_filter(custom,custom_simulated%returns,n=18)
  call assert_true(custom_filtered%ok,'custom model filtering')
  custom%mu_y_parameters=[0.0_dp]
  custom_initial=dnf_filter(custom,custom_simulated%returns,n=12)
  custom_fit=dnf_optimize_custom(custom,custom_simulated%returns,[0.0_dp],custom_setter, &
    n=12,max_iterations=25,tolerance=1.0e-5_dp,calculate_hessian=.false.)
  call assert_true(custom_fit%ok,'custom model optimization')
  call assert_true(custom_fit%log_likelihood>=custom_initial%log_likelihood-1.0e-8_dp, &
    'custom optimizer improves likelihood')

  write(*,'(a)') 'test_models: PASS'
end program test_models
