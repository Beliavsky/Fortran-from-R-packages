! SPDX-License-Identifier: MIT
program asian_monte_carlo
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: requested(4)
  integer :: i

  requested = [character(len=24) :: 'fair_value', 'delta', 'rho', 'vega']
  call bs_malliavin_asian_greeks(100.0_dp, 100.0_dp, 0.03_dp, 1.0_dp, &
    0.25_dp, 0.01_dp, 'call', requested, result, steps=48, paths=20000, seed=7)
  if (result%status /= greeks_ok) error stop trim(result%message)
  do i = 1, size(result%values)
    print '(a24,2x,es16.8,2x,a,es12.4)', trim(result%names(i)), &
      result%values(i), 'MC SE:', result%standard_errors(i)
  end do
end program asian_monte_carlo
