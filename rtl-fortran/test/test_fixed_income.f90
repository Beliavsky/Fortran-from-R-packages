! SPDX-License-Identifier: MIT
program test_fixed_income
  use rtl, only: dp, bond_result, npv_result, irs_result, commodity_weight_result
  use rtl, only: bond_value, npv_value, interest_rate_swap, commodity_swap_prices
  use rtl, only: commodity_swap_from_calendar
  use rtl, only: ymd_to_serial, swap_fut_weight
  implicit none

  type(bond_result) :: bond_output
  type(npv_result) :: npv_output
  type(irs_result) :: irs
  type(commodity_weight_result) :: weight
  real(dp) :: curve_t(5), curve_d(5), swaps(2)
  integer :: trade_date, effective_date, maturity_date, month_start, expiry_date

  bond_output = bond_value(0.05_dp, 0.05_dp, 1.0_dp, 2)
  call assert_true(bond_output%status%ok)
  call assert_close(bond_output%price, 100.0_dp, 1.0e-12_dp)
  call assert_close(bond_output%duration, 0.987804878048780_dp, 1.0e-12_dp)

  curve_t = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
  curve_d = exp(-0.04_dp * curve_t)
  npv_output = npv_value(-375.0_dp, 50.0_dp, 0.5_dp, 250.0_dp, 2.0_dp, &
    curve_t, curve_d, terminal_replaces_cash_flow=.true.)
  call assert_close(npv_output%value, sum(npv_output%present_value), 0.0_dp)
  call assert_close(npv_output%cash_flow(4), 250.0_dp, 0.0_dp)
  npv_output = npv_value(-375.0_dp, 50.0_dp, 0.5_dp, 250.0_dp, 2.0_dp, &
    curve_t, curve_d, terminal_replaces_cash_flow=.false.)
  call assert_close(npv_output%cash_flow(4), 300.0_dp, 0.0_dp)

  trade_date = ymd_to_serial(2020, 1, 4)
  effective_date = ymd_to_serial(2020, 1, 6)
  maturity_date = ymd_to_serial(2022, 1, 6)
  irs = interest_rate_swap(trade_date, effective_date, maturity_date, 1000000.0_dp, &
    "Rec", 0.0_dp, curve_t, spread(1.0_dp, 1, 5), 6, curve_t, spread(1.0_dp, 1, 5), 360)
  call assert_true(irs%status%ok)
  call assert_close(irs%present_value, 0.0_dp, 1.0e-12_dp)
  call assert_true(size(irs%payment_dates) == 4)

  month_start = ymd_to_serial(2020, 9, 1)
  expiry_date = ymd_to_serial(2020, 9, 21)
  weight = swap_fut_weight(month_start, expiry_date)
  call assert_true(weight%status%ok)
  call assert_true(weight%days_first + weight%days_second == 22)
  swaps = commodity_swap_prices([40.0_dp, 41.0_dp], [42.0_dp, 43.0_dp], 0.25_dp)
  call assert_close(swaps(1), 41.5_dp, 0.0_dp)
  call assert_close(swaps(2), 42.5_dp, 0.0_dp)
  swaps = commodity_swap_from_calendar([40.0_dp, 41.0_dp], [42.0_dp, 43.0_dp], &
    month_start, expiry_date)
  call assert_close(swaps(1), 40.0_dp * weight%first_weight + &
    42.0_dp * (1.0_dp - weight%first_weight), 1.0e-14_dp)

  print '(a)', 'test_fixed_income: PASS'

contains

  subroutine assert_close(actual, expected, tolerance)
    real(dp), intent(in) :: actual, expected, tolerance
    if (abs(actual - expected) > tolerance) then
      print '(a,3es24.15)', 'mismatch: ', actual, expected, abs(actual - expected)
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition)
    logical, intent(in) :: condition
    if (.not. condition) error stop 1
  end subroutine assert_true

end program test_fixed_income
