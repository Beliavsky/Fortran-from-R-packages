program ltsa_demo
    use iso_fortran_env, only : int64
    use ltsa, only : dp, ltsa_error, forecast_result, exact_likelihood_result, set_ltsa_seed, &
                     tacvf_arma, dl_simulate, dl_residuals, exact_loglikelihood, trench_forecast
    implicit none
    type(ltsa_error) :: error
    type(forecast_result) :: forecast
    type(exact_likelihood_result) :: likelihood
    real(dp), allocatable :: r(:), z(:), residuals(:)
    real(dp) :: phi(1)

    phi = 0.8_dp
    call tacvf_arma(phi,[real(dp)::],104,1.0_dp,r,error)
    if (.not.error%ok()) error stop trim(error%message)
    call set_ltsa_seed(20260804_int64)
    call dl_simulate(100,r,z,error)
    if (.not.error%ok()) error stop trim(error%message)
    call dl_residuals(r(1:100),z,residuals,error)
    if (.not.error%ok()) error stop trim(error%message)
    likelihood = exact_loglikelihood(r(1:100),z)
    forecast = trench_forecast(z,r,0.0_dp,100,5)

    print '(a,f10.5)', 'exact log-likelihood: ', likelihood%log_likelihood
    print '(a,f10.5)', 'estimated innovation variance: ', likelihood%sigma_sq
    print '(a,f10.5)', 'standardized residual RMS: ', sqrt(sum(residuals**2)/real(size(residuals),dp))
    print '(a,5f10.5)', 'five forecasts: ', forecast%forecasts(1,:)
    print '(a,5f10.5)', 'forecast standard deviations: ', forecast%sd_forecasts(1,:)
end program ltsa_demo
