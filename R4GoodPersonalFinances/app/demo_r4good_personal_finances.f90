! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program demo_r4good_personal_finances
  use r4good_personal_finances
  implicit none
  type(household) :: home
  type(household_member) :: member
  type(portfolio_spec) :: portfolio
  type(portfolio_parameters) :: parameters
  type(lifecycle_result) :: scenario
  real(dp), allocatable :: income(:), essential_spending(:)
  integer :: n, status

  print '(a)', 'R4GoodPersonalFinances modern Fortran demonstration'
  print '(a,f10.4)', 'Survival from 65 to 85: ', &
    gompertz_survival_probability(65.0_dp, 85.0_dp, 80.0_dp, 10.0_dp)
  print '(a,f10.4)', 'Four-percent retirement ruin probability: ', &
    retirement_ruin_probability(0.02_dp, 0.20_dp, 65.0_dp, 88.0_dp, 10.0_dp, 0.04_dp)

  member%name = 'primary'
  member%birth_date = date_from_string('1980-02-15')
  member%mode = 88.0_dp
  member%dispersion = 10.0_dp
  call home%add_member(member, status)
  home%configured_lifespan = 10.0_dp
  home%risk_tolerance = 0.5_dp

  call create_default_portfolio(portfolio)
  portfolio%taxable_accounts = [200000.0_dp, 100000.0_dp]
  portfolio%taxadvantaged_accounts = [50000.0_dp, 25000.0_dp]
  call calculate_portfolio_parameters(portfolio, parameters)
  print '(a,f12.2)', 'Initial financial wealth: ', parameters%value
  print '(a,2f10.4)', 'Initial portfolio weights: ', parameters%weights

  n = int(home%configured_lifespan) + 1
  allocate(income(n), essential_spending(n))
  income = 70000.0_dp
  essential_spending = 45000.0_dp
  call simulate_lifecycle(home, portfolio, date_from_string('2026-08-02'), income, essential_spending, &
    scenario, use_random_returns=.false.)
  print '(a,f12.2)', 'First-year optimal discretionary spending: ', scenario%discretionary_spending(1)
  print '(a,2f10.4)', 'First-year optimal asset weights: ', scenario%total_allocation(:,1)
  print '(a,f12.2)', 'Wealth after first-year return: ', scenario%financial_wealth_end(1)
end program demo_r4good_personal_finances
