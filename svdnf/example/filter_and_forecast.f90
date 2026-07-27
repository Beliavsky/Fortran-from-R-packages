! SPDX-License-Identifier: GPL-3.0-only
program filter_and_forecast
  use svdnf
  implicit none
  type(svm_dynamics) :: dynamics
  type(simulation_result) :: simulated
  type(filter_result) :: filtered
  type(percentile_result) :: filtered_median, predicted_median
  type(forecast_result) :: forecast

  dynamics=dynamics_svm('TaylorWithLeverage',phi=0.96_dp,theta=-1.0_dp, &
    sigma=0.18_dp,rho=-0.45_dp)
  simulated=model_sim(dynamics,t=100,init_vol=-1.0_dp,seed=12345)
  filtered=dnf(dynamics,simulated%returns,n=28)
  filtered_median=extract_vol_perc(filtered,p=0.5_dp)
  predicted_median=extract_vol_perc(filtered,p=0.5_dp,pred=.true.)
  forecast=predict_svdnf(filtered,n_ahead=10,n_sim=1000,confidence=0.95_dp,seed=54321)

  write(*,'(a,f14.6)') 'log likelihood: ',filtered%log_likelihood
  write(*,'(a,f12.6)') 'last filtered median: ', &
    filtered_median%values(size(filtered_median%values))
  write(*,'(a,f12.6)') 'last one-step prediction median: ', &
    predicted_median%values(size(predicted_median%values))
  write(*,'(a,10f10.5)') 'return forecast: ',forecast%mean_return
end program filter_and_forecast
