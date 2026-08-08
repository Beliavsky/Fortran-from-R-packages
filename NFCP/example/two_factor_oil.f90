program two_factor_oil
  use nfcp
  implicit none

  type(nfcp_model_t) :: model
  type(nfcp_futures_simulation_result_t) :: simulation
  type(nfcp_option_result_t) :: option
  real(dp), allocatable :: forecast(:,:), maturities(:,:)
  real(dp) :: initial_state(2), forecast_maturities(4)
  integer :: i, status

  call initialize_model(model, n_factors=2, gbm=.true., n_me=1, n_season=1)
  model%mu = 0.03_dp
  model%mu_rn = 0.01_dp
  model%lambda = [0.02_dp, 0.10_dp]
  model%kappa = [0.0_dp, 1.20_dp]
  model%sigma = [0.20_dp, 0.30_dp]
  model%rho = reshape([1.0_dp, -0.20_dp, -0.20_dp, 1.0_dp], [2, 2])
  model%season_cos = [0.04_dp]
  model%season_sin = [-0.02_dp]
  model%measurement_error = [0.01_dp]

  initial_state = [log(50.0_dp), -0.20_dp]
  forecast_maturities = [0.50_dp, 0.75_dp, 1.25_dp, 2.25_dp]

  call futures_price_forecast(model, initial_state, 0.25_dp, &
                              forecast_maturities, forecast, &
                              percentiles=[0.05_dp, 0.95_dp], status=status)
  if (status /= nfcp_ok) error stop 'futures forecast failed'

  write(*,'(a)') 'Two-factor futures forecast at observation time 0.25'
  write(*,'(a)') ' maturity       median          p05          p95'
  do i = 1, size(forecast,1)
    write(*,'(f9.2,3f13.4)') forecast_maturities(i), forecast(i,:)
  end do

  allocate(maturities(13,3))
  do i = 1, size(maturities,1)
    maturities(i,:) = real(i-1,dp)/12.0_dp + [0.25_dp, 0.75_dp, 1.50_dp]
  end do
  call futures_price_simulate(model, initial_state, 1.0_dp/12.0_dp, maturities, &
                              simulation, seed=2026)
  if (simulation%status /= nfcp_ok) error stop 'futures simulation failed'
  write(*,'(/,a,3f12.4)') 'Final simulated futures: ', simulation%futures_prices(13,:)

  call european_option_value(model, initial_state, futures_maturity=1.0_dp, &
                             option_maturity=0.5_dp, strike=50.0_dp, &
                             risk_free_rate=0.03_dp, is_call=.true., result=option)
  if (option%status /= nfcp_ok) error stop 'option valuation failed'
  write(*,'(a,f12.6)') 'Six-month call value: ', option%value
  write(*,'(a,f12.6)') 'Annualized volatility: ', option%annualized_volatility
end program two_factor_oil
