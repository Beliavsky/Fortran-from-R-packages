! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program lifecycle_simulation
  use r4good_personal_finances
  implicit none
  type(household) :: home
  type(household_member) :: member
  type(portfolio_spec) :: portfolio
  type(lifecycle_result) :: result
  real(dp), allocatable :: income(:), spending(:)
  integer :: n, status, i

  member%name = 'member'
  member%birth_date = date_from_string('1980-02-15')
  member%mode = 88.0_dp
  member%dispersion = 10.0_dp
  call home%add_member(member, status)
  home%configured_lifespan = 8.0_dp
  call create_default_portfolio(portfolio)
  portfolio%taxable_accounts = [200000.0_dp, 100000.0_dp]
  portfolio%taxadvantaged_accounts = [50000.0_dp, 25000.0_dp]
  n = int(home%configured_lifespan) + 1
  allocate(income(n), spending(n))
  income = 70000.0_dp
  spending = 45000.0_dp

  call simulate_lifecycle(home, portfolio, date_from_string('2026-08-02'), income, spending, result, &
    use_random_returns=.true., seed=20260802_i8)
  print '(a)', 'year wealth net_worth discretionary_spending risky_weight'
  do i = 1, result%n_periods
    print '(i4,3f15.2,f12.4)', result%timeline%year(i), result%financial_wealth(i), &
      result%net_worth(i), result%discretionary_spending(i), result%total_allocation(1,i)
  end do
end program lifecycle_simulation
