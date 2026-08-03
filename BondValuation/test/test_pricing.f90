! SPDX-License-Identifier: GPL-3.0-only
program test_pricing
  use bondvaluation
  implicit none
  type(bond_terms) :: terms
  type(bond_schedule) :: schedule
  type(bond_value_result) :: price_result, yield_result
  real(dp) :: expected_zero_dirty, zero_power

  terms%issue_date = date_from_ymd(2020, 1, 15)
  terms%maturity_date = date_from_ymd(2030, 1, 15)
  terms%coupon_frequency = 2
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.0_dp
  terms%day_count_convention = dcc_act_act_icma
  call build_bond_schedule(terms, schedule)
  price_result = bond_price(4.0_dp, date_from_ymd(2024, 4, 15), terms, &
    supplied_schedule=schedule)
  call assert_equal_int(price_result%status, 0, "regular price status")
  call assert_close(price_result%dirty_price, 106.33533492789661_dp, 2.0e-11_dp, &
    "regular dirty price")
  call assert_close(price_result%accrued_interest, 1.25_dp, 2.0e-13_dp, &
    "regular accrued interest")
  call assert_close(price_result%clean_price, 105.08533492789661_dp, 2.0e-11_dp, &
    "regular clean price")
  yield_result = bond_yield(price_result%clean_price, date_from_ymd(2024, 4, 15), &
    terms, supplied_schedule=schedule)
  call assert_equal_int(yield_result%status, 0, "regular yield status")
  call assert_close(yield_result%yield_percent, 4.0_dp, 2.0e-10_dp, &
    "regular yield inversion")

  terms = bond_terms()
  terms%issue_date = date_from_ymd(2013, 11, 30)
  terms%maturity_date = date_from_ymd(2021, 4, 21)
  terms%coupon_frequency = 2
  terms%first_coupon_date = date_from_ymd(2015, 2, 28)
  terms%last_coupon_date = date_from_ymd(2020, 2, 29)
  terms%first_interest_accrual_date = terms%issue_date
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 5.25_dp
  terms%day_count_convention = dcc_act_act_icma
  price_result = bond_price(5.0_dp, date_from_ymd(2014, 10, 15), terms)
  call assert_equal_int(price_result%status, 0, "odd price status")
  call assert_close(price_result%accrued_interest, 4.582872928176794_dp, &
    2.0e-11_dp, "odd accrued interest")
  call assert_close(price_result%dirty_price, 105.81073224263731_dp, &
    2.0e-11_dp, "odd dirty price")
  call assert_close(price_result%clean_price, 101.22785931446052_dp, &
    2.0e-11_dp, "odd clean price")
  call assert_close(price_result%modified_duration_periods, 10.534694027524306_dp, &
    2.0e-10_dp, "odd modified duration")
  call assert_close(price_result%macaulay_duration_periods, 10.798061378212413_dp, &
    2.0e-10_dp, "odd Macaulay duration")
  call assert_close(price_result%convexity_periods, 68.40044885201478_dp, &
    2.0e-10_dp, "odd convexity")
  yield_result = bond_yield(price_result%clean_price, date_from_ymd(2014, 10, 15), &
    terms)
  call assert_close(yield_result%yield_percent, 5.0_dp, 2.0e-10_dp, &
    "odd yield inversion")

  terms = bond_terms()
  terms%issue_date = date_from_ymd(2024, 1, 1)
  terms%maturity_date = date_from_ymd(2025, 1, 1)
  terms%coupon_frequency = 0
  terms%redemption_value = 100.0_dp
  terms%coupon_rate_percent = 0.0_dp
  terms%day_count_convention = dcc_act_act_icma
  call build_bond_schedule(terms, schedule)
  zero_power = discount_coordinate(date_from_ymd(2025, 1, 1), &
    dcc_act_act_icma, terms, schedule) - discount_coordinate(&
    date_from_ymd(2024, 7, 1), dcc_act_act_icma, terms, schedule)
  expected_zero_dirty = 100.0_dp / (1.0_dp + 0.05_dp * zero_power)
  price_result = bond_price(5.0_dp, date_from_ymd(2024, 7, 1), terms, .true., &
    supplied_schedule=schedule)
  call assert_close(price_result%dirty_price, expected_zero_dirty, 2.0e-12_dp, &
    "zero-coupon price")
  call assert_close(price_result%clean_price, expected_zero_dirty, 2.0e-12_dp, &
    "zero-coupon clean price")
  yield_result = bond_yield(price_result%clean_price, date_from_ymd(2024, 7, 1), &
    terms, .true., supplied_schedule=schedule)
  call assert_close(yield_result%yield_percent, 5.0_dp, 2.0e-10_dp, &
    "zero-coupon yield inversion")
  print '(a)', "test_pricing: PASS"

contains

  subroutine assert_close(actual, expected, tolerance, label)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    if (abs(actual - expected) > tolerance) then
      write(*, '(a,2(1x,es24.16))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_equal_int(actual, expected, label)
    integer, intent(in) :: actual, expected
    character(len=*), intent(in) :: label
    if (actual /= expected) then
      write(*, '(a,2(1x,i0))') trim(label)//" mismatch:", actual, expected
      error stop 1
    end if
  end subroutine assert_equal_int

end program test_pricing
