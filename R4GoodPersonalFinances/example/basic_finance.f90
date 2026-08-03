! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program basic_finance
  use r4good_personal_finances
  implicit none
  real(dp), allocatable :: pv(:)
  real(dp) :: cashflows(10)
  integer :: status

  cashflows = 100000.0_dp
  call present_value_stream(cashflows, 0.02_dp, pv, status)

  print '(a,f12.2)', 'Present value of ten annual 100,000 payments: ', pv(1)
  print '(a,f10.4)', 'Optimal risky allocation: ', &
    optimal_risky_asset_allocation(0.05_dp, 0.20_dp, 0.00_dp, 2.0_dp)
  print '(a,f10.2)', 'Purchasing power after 30 years at -2% real: ', &
    purchasing_power(10.0_dp, 30.0_dp, -0.02_dp)
end program basic_finance
