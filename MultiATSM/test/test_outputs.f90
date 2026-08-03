program test_outputs
  use multiatsm_kinds, only : dp
  use multiatsm_types, only : response_result, variance_decomposition_result, forecast_result
  use multiatsm_outputs, only : forecast_yields, impulse_responses, generalized_impulse_responses, &
    forecast_error_variance_decomposition, generalized_fevd, expected_short_rate_component, &
    term_premium, forward_rates
  implicit none
  real(dp) :: phi(2, 2), sigma(2, 2), g0(2, 2), b(3, 2), k0(2), a(3), last(2)
  real(dp) :: states(2, 4), rho0(1), rho1(1, 2), observed(2, 4), yields_forwards(3, 2)
  integer :: maturities(2), forward_maturities(3), info
  real(dp), allocatable :: expected(:, :, :), premiums(:, :), forwards(:, :)
  type(response_result) :: irf, girf
  type(variance_decomposition_result) :: fevd, gfevd
  type(forecast_result) :: forecast

  phi = 0.0_dp
  phi(1, 1) = 0.5_dp
  phi(2, 2) = 0.2_dp
  sigma = 0.0_dp
  sigma(1, 1) = 1.0_dp
  sigma(2, 2) = 4.0_dp
  g0 = 0.0_dp
  g0(1, 1) = 1.0_dp
  g0(2, 2) = 1.0_dp
  b = reshape([1.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, 0.5_dp], [3, 2])
  k0 = [0.1_dp, -0.05_dp]
  a = [0.01_dp, 0.02_dp, 0.03_dp]
  last = [0.2_dp, -0.1_dp]

  call forecast_yields(k0, phi, a, b, last, 3, forecast, info)
  call check(info == 0, 'forecast status')
  call check(maxval(abs(forecast%factors(:, 1) - [0.2_dp, -0.07_dp])) < 1.0e-12_dp, &
    'one-step factor forecast')
  call check(maxval(abs(forecast%yields(:, 1) - (a + matmul(b, forecast%factors(:, 1))))) < 1.0e-12_dp, &
    'yield forecast mapping')

  call impulse_responses(phi, sigma, g0, b, 3, irf, info)
  call check(info == 0, 'IRF status')
  call check(abs(irf%factors(1, 1, 1) - 1.0_dp) < 1.0e-12_dp, 'impact response')
  call check(abs(irf%factors(1, 1, 2) - 0.5_dp) < 1.0e-12_dp, 'dynamic response')
  call generalized_impulse_responses(phi, sigma, g0, b, 3, girf, info)
  call check(info == 0, 'GIRF status')
  call check(abs(girf%factors(2, 2, 1) - 2.0_dp) < 1.0e-10_dp, 'generalized scale')

  call forecast_error_variance_decomposition(phi, sigma, g0, b, 4, fevd, info)
  call check(info == 0, 'FEVD status')
  call check(maxval(abs(sum(fevd%factors(:, :, 4), dim=2) - 1.0_dp)) < 1.0e-12_dp, &
    'factor FEVD normalization')
  call generalized_fevd(phi, sigma, g0, b, 4, gfevd, info)
  call check(info == 0, 'GFEVD status')
  call check(maxval(abs(sum(gfevd%yields(:, :, 4), dim=2) - 1.0_dp)) < 1.0e-12_dp, &
    'yield GFEVD normalization')

  states = reshape([0.02_dp, 0.01_dp, 0.03_dp, 0.015_dp, 0.04_dp, 0.02_dp, 0.05_dp, 0.025_dp], [2, 4])
  rho0 = 0.005_dp
  rho1 = reshape([1.0_dp, 0.0_dp], [1, 2])
  maturities = [1, 2]
  call expected_short_rate_component(k0, phi, states, rho0, rho1, maturities, expected, info, floor_zero=.false.)
  call check(info == 0, 'expected component status')
  call check(abs(expected(1, 1, 1) - 0.025_dp) < 1.0e-12_dp, 'one-period expected component')
  observed(1, :) = expected(1, 1, :) + 0.01_dp
  observed(2, :) = expected(1, 2, :) + 0.02_dp
  call term_premium(observed, expected, 1, premiums, info)
  call check(info == 0, 'term premium status')
  call check(maxval(abs(premiums(1, :) - 0.01_dp)) < 1.0e-12_dp, 'term premium maturity 1')

  forward_maturities = [1, 2, 4]
  yields_forwards(:, 1) = [0.02_dp, 0.03_dp, 0.04_dp]
  yields_forwards(:, 2) = [0.01_dp, 0.015_dp, 0.025_dp]
  call forward_rates(yields_forwards, forward_maturities, forwards, info)
  call check(info == 0, 'forward status')
  call check(abs(forwards(1, 1) - 0.04_dp) < 1.0e-12_dp, 'forward value')
  print '(a)', 'test_outputs: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_outputs
