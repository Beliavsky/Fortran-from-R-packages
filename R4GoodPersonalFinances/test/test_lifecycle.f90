! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program test_lifecycle
  use r4good_personal_finances
  implicit none
  type(household) :: home
  type(household_member) :: member
  type(portfolio_spec) :: portfolio
  type(lifecycle_result) :: result
  type(date_type) :: current_date
  real(dp), allocatable :: income(:), spending(:)
  integer :: st, n

  current_date = date_from_string("2020-07-15")
  member%name = "older"
  member%birth_date = date_from_string("1980-02-15")
  member%mode = 80.0_dp
  member%dispersion = 10.0_dp
  member%max_age = 50.0_dp
  call home%add_member(member, st)
  home%configured_lifespan = 5.0_dp
  home%risk_tolerance = 0.5_dp
  home%consumption_impatience_preference = 0.04_dp
  home%smooth_consumption_preference = 1.0_dp
  call create_default_portfolio(portfolio)
  portfolio%taxable_accounts = [200000.0_dp, 100000.0_dp]
  portfolio%taxadvantaged_accounts = [50000.0_dp, 25000.0_dp]
  n = 6
  allocate(income(n), spending(n))
  income = 30000.0_dp
  spending = 24000.0_dp
  call simulate_lifecycle(home, portfolio, current_date, income, spending, result, use_random_returns=.false., &
    optimizer_max_iterations=3000)
  call assert_true(result%status == r4gpf_success, "lifecycle status")
  call assert_true(result%n_periods == n .and. result%n_assets == 2, "lifecycle dimensions")
  call assert_close(result%financial_wealth(1), 375000.0_dp, 1.0e-8_dp, "initial wealth")
  call assert_close(sum(result%total_allocation(:,1)), 1.0_dp, 1.0e-10_dp, "allocation sum")
  call assert_true(all(result%time_value_discount > 0.0_dp), "discount factors")
  call assert_true(maxval(abs(result%returns(:,1) - portfolio%expected_return(1))) <= 1.0e-14_dp, "expected return path")
  call assert_true(abs(result%financial_wealth_end(1) - result%financial_wealth(1)) > 1.0e-8_dp, "wealth evolves")
  print *, "test_lifecycle: PASS"
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      print *, "FAIL: ", trim(message)
      error stop 1
    end if
  end subroutine assert_true
  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual - expected) > tolerance) then
      print *, "FAIL: ", trim(message), actual, expected
      error stop 1
    end if
  end subroutine assert_close
end program test_lifecycle
