! SPDX-License-Identifier: GPL-3.0-only
program svdnf_demo
  use svdnf
  implicit none
  type(svm_dynamics) :: dynamics
  type(simulation_result) :: simulated
  type(filter_result) :: filtered
  type(percentile_result) :: median
  type(forecast_result) :: forecast

  dynamics = dynamics_svm('Heston',mu=0.04_dp,kappa=3.0_dp,theta=0.04_dp, &
    sigma=0.35_dp,rho=-0.65_dp)
  simulated = model_simulate(dynamics,120,initial_volatility=0.04_dp,seed=20260726)
  filtered = dnf_filter(dynamics,simulated%returns,n=30)
  if (.not. filtered%ok) error stop trim(filtered%message)
  median = extract_vol_percentile(filtered,0.5_dp)
  forecast = predict_filter(filtered,n_ahead=5,n_sim=500,seed=2468)

  write(*,'(a,f14.6)') 'log likelihood: ',filtered%log_likelihood
  write(*,'(a,f12.6)') 'final filtered median: ',median%values(size(median%values))
  write(*,'(a,5f12.6)') 'forecast volatility means: ',forecast%mean_volatility
end program svdnf_demo
