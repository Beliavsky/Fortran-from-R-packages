! SPDX-License-Identifier: MIT
program american_and_implied
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: requested(3)
  real(dp) :: volatility
  integer :: status, iterations, i

  requested = [character(len=24) :: 'fair_value', 'delta', 'gamma']
  call binomial_american_greeks(100.0_dp, 105.0_dp, 0.04_dp, 1.5_dp, &
    0.28_dp, 0.01_dp, 'put', requested, result, steps=500)
  if (result%status /= greeks_ok) error stop trim(result%message)
  do i = 1, size(result%values)
    print '(a24,2x,es16.8)', trim(result%names(i)), result%values(i)
  end do

  call bs_implied_volatility(result%values(1), 100.0_dp, 105.0_dp, 0.04_dp, &
    1.5_dp, 0.01_dp, 'put', volatility, status, iterations)
  print '(a,es16.8,a,i0)', 'European-equivalent implied volatility: ', &
    volatility, ' iterations: ', iterations
end program american_and_implied
