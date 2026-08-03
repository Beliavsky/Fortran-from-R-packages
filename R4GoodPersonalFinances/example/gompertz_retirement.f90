! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program gompertz_retirement
  use r4good_personal_finances
  implicit none
  real(dp) :: mode
  integer :: status

  call gompertz_mode_from_life_expectancy(86.0_dp, 25.0_dp, 8.88_dp, mode, status, 115.0_dp)
  print '(a,f10.4)', 'Calibrated Gompertz mode: ', mode
  print '(a,f10.4)', 'Life expectancy: ', life_expectancy(25.0_dp, mode, 8.88_dp, 115.0_dp)
  print '(a,f10.6)', 'Retirement ruin probability: ', &
    retirement_ruin_probability(0.02_dp, 0.20_dp, 65.0_dp, 88.0_dp, 10.0_dp, 0.04_dp)
end program gompertz_retirement
