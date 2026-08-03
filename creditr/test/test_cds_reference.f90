! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program test_cds_reference
  use creditr
  implicit none
  type(rate_quote_t), allocatable :: quotes(:)
  type(conventions_t) :: conventions
  type(zero_curve_t) :: curve
  type(cds_contract_t) :: contract
  type(cds_result_t) :: result
  integer :: status, failures, i
  real(kind=dp), parameter :: spreads(4) = [12354.529_dp, 243.28_dp, 9106.8084_dp, 1737.7289_dp]
  real(kind=dp), parameter :: coupons(4) = [500.0_dp, 100.0_dp, 500.0_dp, 500.0_dp]
  real(kind=dp), parameter :: dirty_ref(4) = [5704737.09949105_dp, 644893.448654649_dp, &
    5606643.24973030_dp, 3216555.65224484_dp]
  real(kind=dp), parameter :: clean_ref(4) = [5742237.09949105_dp, 652393.448654649_dp, &
    5644143.24973030_dp, 3254055.65224484_dp]
  real(kind=dp), parameter :: hazard_ref(4) = [2.09219959185728_dp, 0.041005885054131344_dp, &
    1.540780446_dp, 0.29314536913622602_dp]
  real(kind=dp), parameter :: spread_dv01_ref(4) = [21.7465214447_dp, 4282.0493381156_dp, &
    42.2023323154_dp, 1574.8583824976_dp]
  real(kind=dp), parameter :: recovery_ref(4) = [-95273.6481793093_dp, -1109.6459348547_dp, &
    -93158.5959274128_dp, -30806.8298468545_dp]

  failures = 0
  call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
  call add_conventions('USD', conventions, status)
  call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)

  contract%trade_date = make_date(2014, 4, 15)
  contract%maturity = make_date(2019, 6, 20)
  contract%use_maturity = .true.
  contract%currency = 'USD'
  contract%recovery = 0.4_dp
  contract%notional = 1.0e7_dp

  do i = 1, 4
    contract%spread_bps = spreads(i)
    contract%coupon_bps = coupons(i)
    result = price_cds(contract, curve, status, quotes)
    call check(status == creditr_ok, 'CDS status')
    call check_relative(result%upfront, dirty_ref(i), 2.0e-6_dp, 'dirty upfront')
    call check_relative(result%principal, clean_ref(i), 2.0e-6_dp, 'clean upfront')
    call check_relative(result%spread_dv01, spread_dv01_ref(i), 2.0e-5_dp, 'spread DV01')
    call check_relative(result%recovery_risk_01, recovery_ref(i), 2.0e-5_dp, 'recovery risk 01')
    if (i /= 3) call check_relative(result%hazard_rate, hazard_ref(i), 1.0e-5_dp, 'hazard rate')
    call check_close(result%accrual, -37500.0_dp * coupons(i) / 500.0_dp, 1.0e-5_dp, 'accrual')
  end do

  if (failures /= 0) error stop 'test_cds_reference failed'
  print '(a)', 'test_cds_reference: PASS'

contains

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      failures = failures + 1
      print '(a)', 'FAIL: ' // trim(label)
    end if
  end subroutine check

  subroutine check_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call check(abs(actual - expected) <= tolerance, label)
  end subroutine check_close

  subroutine check_relative(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call check(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)), label)
  end subroutine check_relative

end program test_cds_reference
