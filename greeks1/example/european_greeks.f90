! SPDX-License-Identifier: MIT
program european_greeks
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: requested(6)
  integer :: i

  requested = [character(len=24) :: 'fair_value', 'delta', 'gamma', &
    'vega', 'theta', 'rho']
  call bs_european_greeks(120.0_dp, 100.0_dp, 0.02_dp, 4.5_dp, 0.22_dp, &
    0.015_dp, 'put', requested, result)
  if (result%status /= greeks_ok) error stop trim(result%message)
  do i = 1, size(result%values)
    print '(a24,2x,es16.8)', trim(result%names(i)), result%values(i)
  end do
end program european_greeks
