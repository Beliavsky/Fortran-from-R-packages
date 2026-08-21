program demo_tsa
  use tsa
  use tseries_random, only : seed_random
  implicit none

  real(dp), allocatable :: x(:)
  type(arimax_result) :: fit
  type(tsa_test_result) :: lb

  call seed_random(12345)
  call arima_sim([0.65_dp], [0.25_dp], 0, 500, x, sigma=0.4_dp)
  fit = arima_fit(x, 1, 0, 1)
  lb = lb_test(fit%residuals, 12, 2)

  write(*,'(a,f9.4)') 'estimated AR(1): ', fit%ar(1)
  write(*,'(a,f9.4)') 'estimated MA(1): ', fit%ma(1)
  write(*,'(a,f9.4)') 'innovation variance: ', fit%sigma2
  write(*,'(a,f9.4)') 'Ljung-Box p-value: ', lb%p_value
end program demo_tsa
