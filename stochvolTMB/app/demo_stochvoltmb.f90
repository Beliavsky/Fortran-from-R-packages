program demo_stochvoltmb
  use stochvoltmb
  implicit none
  type(sv_rng_state) :: rng
  type(sv_parameters) :: true_params
  type(sv_simulation) :: sim
  type(sv_fit_result) :: fit
  type(sv_control) :: control
  type(sv_prediction) :: prediction
  type(sv_prediction_summary) :: forecast
  real(dp) :: probabilities(2)

  true_params%sigma_y=0.30_dp
  true_params%sigma_h=0.25_dp
  true_params%phi=0.92_dp
  call rng%seed(20260804_8)
  call sim_sv(true_params,160,sv_gaussian,rng,sim)

  control%max_outer_iter=140
  control%compute_covariance=.true.
  call estimate_parameters(sim%y,sv_gaussian,fit,control=control)

  write(*,'(a,l1)') 'converged: ',fit%converged
  write(*,'(a,f10.5)') 'sigma_y:  ',fit%params%sigma_y
  write(*,'(a,f10.5)') 'sigma_h:  ',fit%params%sigma_h
  write(*,'(a,f10.5)') 'phi:      ',fit%params%phi
  write(*,'(a,f12.4)') 'logLik:   ',fit%log_likelihood

  call predict_sv(fit,5,2000,rng,prediction,include_parameters=.true.)
  probabilities=[0.025_dp,0.975_dp]
  call summarize_prediction(prediction,probabilities,forecast,include_mean=.true.)
  write(*,'(/,a)') 'Five-step predictive return summary:'
  write(*,'(a)') ' step        mean        q025        q975'
  call print_forecast(forecast)
contains
  subroutine print_forecast(s)
    type(sv_prediction_summary), intent(in) :: s
    integer :: i
    do i=1,size(s%y_mean)
      write(*,'(i5,3f12.6)') i,s%y_mean(i),s%y_quantiles(i,1),s%y_quantiles(i,2)
    end do
  end subroutine print_forecast
end program demo_stochvoltmb
