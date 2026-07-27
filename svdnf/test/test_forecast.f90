! SPDX-License-Identifier: GPL-3.0-only
program test_forecast
  use svdnf
  use test_support
  implicit none
  type(svm_dynamics) :: dynamics
  type(simulation_result) :: simulated
  type(filter_result) :: filtered
  type(forecast_result) :: forecast1, forecast2

  dynamics=dynamics_svm('TaylorWithLeverage',phi=0.94_dp,theta=-1.1_dp,sigma=0.18_dp,rho=-0.45_dp)
  simulated=model_simulate(dynamics,70,initial_volatility=-1.1_dp,seed=777)
  filtered=dnf_filter(dynamics,simulated%returns,n=22)
  call assert_true(filtered%ok,'filter before forecast')
  forecast1=predict_filter(filtered,n_ahead=8,n_sim=300,confidence=0.9_dp,seed=1357)
  forecast2=predict_filter(filtered,n_ahead=8,n_sim=300,confidence=0.9_dp,seed=1357)
  call assert_true(forecast1%ok .and. forecast2%ok,'forecast status')
  call assert_vector_close(forecast1%mean_volatility,forecast2%mean_volatility,0.0_dp, &
    'deterministic seeded forecast')
  call assert_true(all(forecast1%lower_volatility<=forecast1%mean_volatility),'volatility lower interval')
  call assert_true(all(forecast1%upper_volatility>=forecast1%mean_volatility),'volatility upper interval')
  call assert_true(all(forecast1%lower_return<=forecast1%mean_return),'return lower interval')
  call assert_true(all(forecast1%upper_return>=forecast1%mean_return),'return upper interval')

  write(*,'(a)') 'test_forecast: PASS'
end program test_forecast
