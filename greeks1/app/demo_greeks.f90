! SPDX-License-Identifier: MIT
program demo_greeks
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: requested(5)
  integer :: i

  requested = [character(len=24) :: 'fair_value', 'delta', 'gamma', &
    'vega', 'theta']
  call option_greeks(100.0_dp, 100.0_dp, 0.02_dp, 1.0_dp, 0.30_dp, &
    0.01_dp, 'Black_Scholes', 'European', 'call', requested, result)
  if (result%status /= greeks_ok) error stop trim(result%message)
  print '(a)', 'European Black-Scholes call'
  do i = 1, size(result%values)
    print '(2x,a20,2x,es16.8)', trim(result%names(i)), result%values(i)
  end do
end program demo_greeks
